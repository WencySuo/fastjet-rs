use log::warn;
use std::panic;

use crate::constants::PI;
use crate::constants::TWO_PI;

#[cfg(not(feature = "simd"))]
use crate::proxy_jet::BriefJet;
use crate::proxy_jet::EEBriefJet;
use crate::proxy_jet::ProxyJet;
use crate::proxy_jet::Tile;
use crate::proxy_jet::TiledJet;
use crate::pseudo_jet::PseudoJet;
#[cfg(feature = "simd")]
use std::simd::prelude::*;
#[cfg(feature = "simd")]
use std::simd::{Select, Simd, StdFloat};

// TODO: add automatic strategies for clustering like BEST: N2Plain, N2Tiling
pub struct ClusterSequence {
    pub particles: Vec<PseudoJet>,
    pub jetdef: JetDefinition,
    _invr2: f64,
    r2: f64,
    _history: Vec<HistoryElement>,
    init_n: usize,
    //TODO: find if we can move tiled jets related vars out of here
    tiles_struct: ClusterSequenceTiles,
}

pub struct ClusterBuffer {
    ptr: *mut f64,
    n: usize,
}

impl ClusterBuffer {
    pub fn new(n: usize) -> Self {
        // We need 6 arrays: rap, phi, kt2, nn_dist, di_j, nn_idx
        // We must ensure 'n' is padded to be a multiple of LANES
        // so that the NEXT array also starts on a 128-byte boundary.

        // Total size: (padded_n * 6 arrays * 8 bytes per element)
        // Alignment: 128 bytes (0x80)
        let layout = std::alloc::Layout::from_size_align(n * 6 * 8, 64).unwrap();
        let ptr = unsafe { std::alloc::alloc_zeroed(layout) as *mut f64 };

        if ptr.is_null() {
            panic!("Allocation failed");
        }

        Self { ptr, n }
    }

    pub fn get_context(&mut self) -> ClusterContext<'_> {
        let p_n = self.n; // Use the padded size for offsets!
        unsafe {
            ClusterContext {
                rap: std::slice::from_raw_parts_mut(self.ptr, self.n),
                phi: std::slice::from_raw_parts_mut(self.ptr.add(p_n), self.n),
                kt2: std::slice::from_raw_parts_mut(self.ptr.add(p_n * 2), self.n),
                nn_dist: std::slice::from_raw_parts_mut(self.ptr.add(p_n * 3), self.n),
                di_j: std::slice::from_raw_parts_mut(self.ptr.add(p_n * 4), self.n),
                nn_idx: std::slice::from_raw_parts_mut(self.ptr.add(p_n * 5) as *mut u64, self.n),
            }
        }
    }
}

impl Drop for ClusterBuffer {
    fn drop(&mut self) {
        let layout = std::alloc::Layout::from_size_align(self.n * 6 * 8, 64).unwrap();
        unsafe {
            std::alloc::dealloc(self.ptr as *mut u8, layout);
        }
    }
}

pub struct ClusterContext<'a> {
    pub rap: &'a mut [f64], // Changed to mut
    pub phi: &'a mut [f64], // Changed to mut
    pub kt2: &'a mut [f64], // Changed to mut
    pub nn_dist: &'a mut [f64],
    pub di_j: &'a mut [f64],
    pub nn_idx: &'a mut [u64],
}

#[derive(Default)]
pub struct ClusterSequenceTiles {
    _tiles: Vec<Tile>,
    _tiles_eta_min: f64,
    _tiles_eta_max: f64,
    _tile_size_eta: f64,
    _tile_size_phi: f64,
    _n_tiles_phi: isize,
    _tiles_ieta_min: isize,
    _tiles_ieta_max: isize,
}

pub struct JetDefinition {
    pub algorithm: Algorithm,
    pub r: f64,
    pub scheme: RecombinationScheme,
    pub strategy: Strategy,
    pub extra_param: Option<f64>,
}

pub enum Strategy {
    Best,
    N2Plain,
    N2Tiling,
    N2PlainEEAccurate,
    N2PlainEE,
}

#[derive(PartialEq, Eq)]
pub enum Algorithm {
    AntiKt,
    Kt,
    Cambridge,
    GenKt,
    EeKt,
    EeGenKt,
}

pub enum RecombinationScheme {
    EScheme,
    PTScheme,
}

use RecombinationScheme::*;

impl JetDefinition {
    pub fn new(
        algorithm: Algorithm,
        r: f64,
        scheme: RecombinationScheme,
        strategy: Strategy,
        extra_param: Option<f64>,
    ) -> Self {
        JetDefinition {
            algorithm,
            r,
            scheme,
            strategy,
            extra_param,
        }
    }

    //TODO: should add reset function to make it faster
    pub fn recombine(&self, jet_a: &PseudoJet, jet_b: &PseudoJet) -> PseudoJet {
        match self.scheme {
            EScheme => jet_b + jet_a,
            _ => {
                panic!("scheme not yet supported")
            }
        }
    }

    //TODO: acc add different preprocess for Other recombine algorithms
    pub fn preprocess(&self, _jet_a: &mut PseudoJet) {
        //E-scheme does not need any preprocessing
        match self.scheme {
            EScheme => {}
            _ => {
                panic!("scheme not yet supported")
            }
        }
    }
}

// -----------ENUMS--------------

pub const INVALID: usize = usize::MAX;
pub const BEAMJET: usize = usize::MAX - 1;
pub const INEXISTENT_PARENT: usize = usize::MAX - 2;

pub struct HistoryElement {
    parent1: usize,
    parent2: usize,
    child: usize,
    jet_index: usize,
    d_ij: f64,
    curr_max_d_ij: f64,
}

impl ClusterSequence {
    // -----------CONSTRUCTOR--------------
    //TODO: verify if this is right place for constant
    pub fn new(mut particles: Vec<PseudoJet>, jetdef: JetDefinition) -> Self {
        let r2 = jetdef.r * jetdef.r;
        let _invr2 = 1.0 / r2;
        //recomb will produce 2x particles len in max case
        let _history: Vec<HistoryElement> = Vec::with_capacity(particles.len() * 2);
        let init_n = particles.len();
        particles.reserve(particles.len());
        // filling init history requires
        let mut clust_seq = ClusterSequence {
            particles,
            jetdef,
            _invr2,
            r2,
            _history,
            init_n,
            tiles_struct: ClusterSequenceTiles::default(),
        };

        clust_seq.fill_init_history();

        clust_seq
    }

    #[inline]
    fn jet_scale_for_algorithm(
        jet: &PseudoJet,
        algorithm: &Algorithm,
        extra_param: Option<f64>,
    ) -> f64 {
        let mut kt2 = jet.kt2();
        match algorithm {
            Algorithm::Kt => kt2,
            Algorithm::AntiKt => {
                if kt2 > 1e-300 {
                    1. / kt2
                } else {
                    1e300
                }
            }
            Algorithm::Cambridge => 1.0,
            Algorithm::GenKt => {
                let p = extra_param.unwrap_or(0.0);
                // TODO: figure out if there's a better way to handle this
                if p <= 0.0 && kt2 < 1e-300 {
                    kt2 = 1e-300;
                }
                kt2.powf(p)
            }
            Algorithm::EeKt => jet.e() * jet.e(),
            Algorithm::EeGenKt => {
                let p = extra_param.unwrap();
                let e = jet.e();
                let mut scale = e * e;
                if p <= 0.0 && scale < 1e-300 {
                    scale = 1e-300;
                    kt2 = scale.powf(p);
                }
                kt2
            }
        }
    }

    #[inline]
    fn n_particles(&self) -> usize {
        self.init_n
    }
    //TODO: check if _decant (which is transfering jet_def and other info)
    // into clust_seq internal variables is acc necessary? some pointer wizardry goin on there
    pub fn initialize_and_run_no_decant(&mut self) {
        //event is empty so exit
        if self.n_particles() == 0 {
            return;
        }

        if self.jetdef.algorithm == Algorithm::EeGenKt || self.jetdef.algorithm == Algorithm::EeKt {
            if self.jetdef.algorithm == Algorithm::EeKt {
                assert!(self.jetdef.r > PI);
                self._invr2 = 1.0;
            } else {
                if self.jetdef.r > PI {
                    self.r2 = 2.0 * (3.0 + (self.jetdef.r).cos());
                } else {
                    self.r2 = 2.0 * (1.0 - (self.jetdef.r).cos());
                }
                self._invr2 = 1.0 / self.r2;
            }
        }

        //TODO: when implementing other strategies have to handle here
        // since we only implemented kt and N2 this func is pretty naive

        match self.jetdef.strategy {
            Strategy::Best => {
                //TODO: implement best strategy
            }
            Strategy::N2Plain => {
                //TODO: simd can only be used with BriefJets that are not EE
                // explicitly include these checks in the future
                #[cfg(feature = "simd")]
                self.simple_n2_cluster_simd();

                #[cfg(not(feature = "simd"))]
                self.simple_n2_cluster::<BriefJet>();
            }
            Strategy::N2Tiling => {
                self.tiled_n2_cluster();
            }
            Strategy::N2PlainEEAccurate => {
                // TODO: decide how to pass down flag for using accurate bj dist
                self.simple_n2_cluster::<EEBriefJet>();
            }
            Strategy::N2PlainEE => {
                self.simple_n2_cluster::<EEBriefJet>();
            }
        }
    }

    fn fill_init_history(&mut self) {
        //get total energy of the event
        let mut tot_e = 0.0;

        self.particles.iter_mut().enumerate().for_each(|(i, jet)| {
            self._history.push(HistoryElement {
                parent1: INEXISTENT_PARENT,
                parent2: INEXISTENT_PARENT,
                child: INVALID,
                jet_index: i,
                d_ij: 0.0,
                curr_max_d_ij: 0.0,
            });

            self.jetdef.preprocess(jet);
            jet.set_cluster_hist_index(i);

            tot_e += jet.e();
        });
        self.init_n = self.particles.len();
    }

    // -----------GETTERS--------------
    pub fn history(&self) -> &Vec<HistoryElement> {
        &self._history
    }

    pub fn inclusive_jets(&self, pmin: f64) -> Vec<&PseudoJet> {
        let dcut = pmin * pmin;
        match self.jetdef.algorithm {
            Algorithm::Kt => {
                // loop through history in reverse
                self._history
                    .iter()
                    .rev()
                    .take_while(|hist_elem| hist_elem.curr_max_d_ij >= dcut)
                    .filter_map(|hist_elem| {
                        if hist_elem.parent2 == BEAMJET && hist_elem.d_ij >= dcut {
                            Some(&self.particles[self._history[hist_elem.parent1].jet_index])
                        } else {
                            None
                        }
                    })
                    .collect()
            }
            Algorithm::EeGenKt | Algorithm::EeKt | Algorithm::GenKt | Algorithm::AntiKt => self
                ._history
                .iter()
                .rev()
                .filter_map(|hist_elem| {
                    if hist_elem.parent2 == BEAMJET {
                        let parent1 = hist_elem.parent1;
                        let jet = &self.particles[self._history[parent1].jet_index];
                        if jet.perp2() >= dcut { Some(jet) } else { None }
                    } else {
                        None
                    }
                })
                .collect(),
            Algorithm::Cambridge => self
                ._history
                .iter()
                .rev()
                .take_while(|hist_elem| hist_elem.parent2 == BEAMJET)
                .filter_map(|hist_elem| {
                    let parent1 = hist_elem.parent1;
                    let jet = &self.particles[self._history[parent1].jet_index];
                    if jet.perp2() >= dcut { Some(jet) } else { None }
                })
                .collect(),
        }
    }

    // note: exclusive jets are not usable well with antikt (only testable later)
    pub fn get_n_exclusive_jets(&self, dcut: f64) -> usize {
        // max_dij so far tells us where in the loop a jet < dcut
        let jets_above_dcut = self
            ._history
            .iter()
            .rev()
            .take_while(|hist_elem| hist_elem.curr_max_d_ij > dcut)
            .count();
        // now we return total jets - jets still in list (particles.len() - jets_above_dcut)
        2 * self.init_n - (self.particles.len() - jets_above_dcut)
    }

    pub fn exclusive_jets_dcut(&self, dcut: f64) -> Vec<&PseudoJet> {
        self.exclusive_jets(self.get_n_exclusive_jets(dcut))
    }

    pub fn exclusive_jets(&self, n_jets: usize) -> Vec<&PseudoJet> {
        //check if we are requesting more particles than possible
        if n_jets > self.init_n {
            panic!(
                "Requesting {} jets which are more than the {} particles in the event",
                n_jets, self.init_n
            )
        }
        self.exclusive_jets_up_to(n_jets)
    }

    pub fn exclusive_jets_up_to(&self, n_jets: usize) -> Vec<&PseudoJet> {
        // warnings for exclusive jets using algorithms that don't support it explicitly.
        let algorithm = &self.jetdef.algorithm;
        if *algorithm != Algorithm::Cambridge
            && *algorithm != Algorithm::Kt
            && *algorithm != Algorithm::GenKt
            && *algorithm != Algorithm::EeKt
            && !((*algorithm == Algorithm::GenKt) || *algorithm == Algorithm::EeGenKt)
            && (self.jetdef.extra_param.unwrap_or(0.0) >= 0.0)
        {
            warn!(
                "dcut and exclusive jets for jet-finders other than kt, C/A or genkt with p>=0 should be interpreted with care."
            );
        }

        // recombining jets means we lose on jet each merge
        // idx in histortical clusterings to stop on
        let mut stop_idx = 2 * self.init_n - n_jets;

        //check if we are requesting more particles than possible

        // max exclusive jets are constrained by # of init input partilces
        if stop_idx < self.init_n {
            stop_idx = self.init_n
        }

        let mut out = Vec::new();

        // loop thru hist to check if parent is last particle before stop_idx
        // also jet types are near usize::max so they will always > stop_idx
        for h in self._history.iter().skip(stop_idx) {
            if h.parent1 < stop_idx {
                out.push(&self.particles[self._history[h.parent1].jet_index]);
            }
            if h.parent2 < stop_idx {
                out.push(&self.particles[self._history[h.parent2].jet_index]);
            }
        }

        out
    }

    //TODO: in fastjet this is a seperate class so should probably move also
    pub fn min_max_rap(&self) -> (f64, f64, f64) {
        //TODO why are these magic numbers here???
        // seems like its the max nrap binning we want to consider?
        // in this config each bin has width of 1 rap
        let nrap = 20;
        let nbins = nrap * 2;
        let mut min_rap = f64::MAX;
        let mut max_rap = -f64::MAX;
        let mut counts = vec![0; nbins];

        let mut ubin: usize;
        for particle in &self.particles {
            //when rap is infinity then ignore the particle
            if particle.E() == particle.pz().abs() {
                continue;
            }
            let rap = particle.rap();
            if rap < min_rap {
                min_rap = rap;
            }
            if rap > max_rap {
                max_rap = rap;
            }

            // rap will be a float so cast to i64 for binning
            let ibin = (rap + nrap as f64) as i64;
            //the leftmost and rightmost bins will be used for overflow
            // thus bin =0 goes from [-infty, -19] and bin=39 [20, infty]
            if ibin < 0 {
                ubin = 0;
            } else if ibin >= nbins as i64 {
                ubin = nbins - 1;
            } else {
                ubin = ibin as usize;
            }
            counts[ubin] += 1;
        }

        //given all binning find count of busiest bin
        // particles will never be empty so just unwrap
        let busiest_bin = counts.iter().max().unwrap();
        //more magic numbers
        let allowed_max_fraction: f64 = 0.25;
        let min_multiplicity: f64 = 4.0;

        //find how much we can stuff in edge bins
        let mut allowed_max_cumul = min_multiplicity
            .max(*busiest_bin as f64 * allowed_max_fraction)
            .floor();
        // chance min_mult is more than busiest bin count
        // in which case allowed_max_cumul should not be greater than busiest bin count
        if allowed_max_cumul > *busiest_bin as f64 {
            allowed_max_cumul = *busiest_bin as f64;
        }

        // scan rap bins from left to find min rap for our tiling
        let mut cum_lo = 0.0;
        let mut cumul_2 = 0.0; //some internal variable class
        ubin = 0;
        for (i, &cnt) in counts.iter().enumerate().take(nbins) {
            cum_lo += cnt as f64;
            if cum_lo >= allowed_max_cumul {
                let y = (i as f64) - (nrap as f64);
                if y > min_rap {
                    min_rap = y;
                    ubin = i;
                }
                break;
            }
        }
        //TODO prolly add an assert that we actally found a bin
        assert!(ubin != nbins);
        cumul_2 += cum_lo * cum_lo;
        let ibin_lo = ubin;

        let mut cum_hi = 0.0;
        // do same as right except from hi to lo
        ubin = 0;
        for i in (0..nbins).rev() {
            cum_hi += counts[i] as f64;
            if cum_hi >= allowed_max_cumul {
                //RHS side of last bin
                let y = i as f64 - nrap as f64 + 1.0;
                if y < max_rap {
                    max_rap = y;
                }
                ubin = i;
                break;
            }
        }

        assert!(ubin != 0);

        let ibin_hi = ubin;

        assert!(ibin_hi >= ibin_lo);

        //cumul2 is the sum of all contents
        if ibin_lo == ibin_hi {
            cumul_2 = (cum_lo + cum_hi - counts[ibin_lo] as f64).powi(2);
        } else {
            cumul_2 += cum_hi * cum_hi;
            //find rest of bins
            for &cnt in counts.iter().take(ibin_hi).skip(ibin_lo + 1) {
                cumul_2 += cnt as f64 * cnt as f64;
            }
        }

        (min_rap, max_rap, cumul_2)
    }

    /// The neighbourhood of a tile is set up as follows:
    /// ```text
    ///       LRR
    ///       LXR
    ///       LLR
    ///
    /// tiles contains XLLLLRRRR with pointers
    ///                     |   \ RH_tiles
    ///                     \ surrounding_tiles
    /// ```
    //setup correct tile size and point tiles to their neighbors
    pub fn initialize_tiles(&mut self) {
        // want to bound this for very low DeltaR to avoid memory blow ups
        let default_tile_size: f64 = self.jetdef.r.max(0.1);
        self.tiles_struct._tile_size_eta = default_tile_size;

        // for phi we need to check spacing against 2pi
        // when phi is <3 no tiling is done so min = 3
        self.tiles_struct._n_tiles_phi = ((TWO_PI / default_tile_size).floor() as isize).max(3);
        self.tiles_struct._tile_size_phi = (TWO_PI) / self.tiles_struct._n_tiles_phi as f64;

        // find the min and max rap/eta that we should be using for this analysis
        let (min_rap, max_rap, _cumul_2) = self.min_max_rap();

        //find min/max values for these tiles and figure out why _tiles_eta_min isnt just min_rap?
        self.tiles_struct._tiles_ieta_min =
            (min_rap / self.tiles_struct._tile_size_eta).floor() as isize;
        self.tiles_struct._tiles_ieta_max =
            (max_rap / self.tiles_struct._tile_size_eta).floor() as isize;
        self.tiles_struct._tiles_eta_min =
            self.tiles_struct._tiles_ieta_min as f64 * self.tiles_struct._tile_size_eta;
        self.tiles_struct._tiles_eta_max =
            self.tiles_struct._tiles_ieta_max as f64 * self.tiles_struct._tile_size_eta;

        // allocate the vector for the size of tiles before pushing
        // TODO: figure out how to access internal _tiles struct with object
        // (should it just be defined as TiledJet object and no trait for PseudoJet)
        // TODO: C++ code does not need to init vectors to get pointers could be perf loss
        let mut tiles: Vec<Tile> = Vec::new();
        tiles.resize(
            (self.tiles_struct._tiles_ieta_max - self.tiles_struct._tiles_ieta_min + 1) as usize
                * self.tiles_struct._n_tiles_phi as usize,
            Tile::new(),
        );

        // now link all these tiles together
        // DEBUGGING TODO: check that the tile_idx logic makes sense and we actually don't get empty or overwritten begin_tiles[]
        for ieta in self.tiles_struct._tiles_ieta_min..=self.tiles_struct._tiles_ieta_max {
            for iphi in 0..self.tiles_struct._n_tiles_phi {
                let mut tile_idx = 0;
                let tile: &mut Tile = &mut tiles[self._tile_int_index(ieta, iphi)];

                // first tile has no HEAD
                tile.head = Option::None;

                //begin_tiles ppoints to neighboring tiles including itself?
                tile.begin_tiles[tile_idx] = self._tile_int_index(ieta, iphi);
                tile_idx += 1;

                // after first loop of ieta we can now check prev column
                // now left right middle
                if ieta > self.tiles_struct._tiles_ieta_min {
                    for idphi in [-1isize, 0, 1] {
                        tile.begin_tiles[tile_idx] = self._tile_int_index(ieta - 1, iphi + idphi);
                        tile_idx += 1;
                    }
                }

                // note: phi doesnt need bound checks since it mod 2pi
                // LHS includes under curr elem
                tile.begin_tiles[tile_idx] = self._tile_int_index(ieta, iphi - 1);
                tile_idx += 1;

                // first RHS is above curr elem
                let rh_tile_start = tile_idx;
                tile.begin_tiles[tile_idx] = self._tile_int_index(ieta, iphi + 1);
                tile_idx += 1;

                // only set last R if we are not at max
                if ieta < self.tiles_struct._tiles_ieta_max {
                    for idphi in [-1isize, 0, 1] {
                        tile.begin_tiles[tile_idx] = self._tile_int_index(ieta + 1, iphi + idphi);
                        tile_idx += 1;
                    }
                }

                // ieta > self._tiles_ieta_min  and when ieta < self._tiles_ieta_max
                tile.begin_len = tile_idx;
                tile.rh_tiles = rh_tile_start..tile_idx;
            }
        }

        // assign tiles to self.tiles_struct._tiles
        self.tiles_struct._tiles = tiles;
    }

    // get index even if iphi is negative since wraparound
    fn _tile_int_index(&self, ieta: isize, iphi: isize) -> usize {
        // use mod to get wraparound behavior but -1 mod n is == -1 so add by n
        ((ieta - self.tiles_struct._tiles_ieta_min) * self.tiles_struct._n_tiles_phi
            + ((iphi + self.tiles_struct._n_tiles_phi) % self.tiles_struct._n_tiles_phi))
            as usize
    }

    fn _tile_index(&self, eta: f64, phi: f64) -> usize {
        let mut ieta: isize;

        //bound eta into an int
        if eta <= self.tiles_struct._tiles_eta_min {
            ieta = 0;
        } else if eta >= self.tiles_struct._tiles_eta_max {
            // this diff is the full range of ieta why not just use ntiles eta???
            // idk you tell me
            ieta = self.tiles_struct._tiles_ieta_max - self.tiles_struct._tiles_ieta_min;
        } else {
            ieta = ((eta - self.tiles_struct._tiles_eta_min) / self.tiles_struct._tile_size_eta)
                as isize;
            //protect against casting
            if ieta > self.tiles_struct._tiles_ieta_max - self.tiles_struct._tiles_ieta_min {
                ieta = self.tiles_struct._tiles_ieta_max - self.tiles_struct._tiles_ieta_min
            }
        };

        let iphi: isize = ((phi + 2.0 * PI) / self.tiles_struct._tile_size_phi) as isize
            % self.tiles_struct._n_tiles_phi;
        (iphi + ieta * self.tiles_struct._n_tiles_phi) as usize
    }

    pub fn _tj_set_jetinfo(
        &mut self,
        jets: &mut Vec<TiledJet>,
        tiled_jet: Option<TiledJet>,
        particle_index: usize,
    ) -> usize {
        let mut jet = if let Some(tj) = tiled_jet {
            tj
        } else {
            let pseudo_jet = &self.particles[particle_index];
            TiledJet::_bj_set_jetinfo(
                particle_index,
                pseudo_jet,
                self.r2,
                ClusterSequence::jet_scale_for_algorithm(
                    pseudo_jet,
                    &self.jetdef.algorithm,
                    self.jetdef.extra_param,
                ),
                particle_index,
            )
        };

        let tile_index = self._tile_index(jet.eta(), jet.phi());
        jet.tile_index = tile_index;

        let old_head = self.tiles_struct._tiles[tile_index].head;
        jet.prev_jet = None;
        jet.next_jet = old_head;

        let new_jet_id = jets.len();
        jets.push(jet);

        if let Some(next) = old_head {
            jets[next].prev_jet = Some(new_jet_id);
        }

        self.tiles_struct._tiles[tile_index].head = Some(new_jet_id);
        new_jet_id
    }

    pub fn bj_remove_from_tiles(&mut self, jets: &mut [TiledJet], jet_idx: usize) {
        let tile = &mut self.tiles_struct._tiles[jets[jet_idx].tile_index];

        // if prev is null then jet is head
        if let Some(prev) = jets[jet_idx].prev_jet {
            // link prev pointer to next pointer to remove jet
            jets[prev].next_jet = jets[jet_idx].next_jet;
        } else {
            tile.head = jets[jet_idx].next_jet;
        }

        if let Some(next) = jets[jet_idx].next_jet {
            jets[next].prev_jet = jets[jet_idx].prev_jet;
        }
    }

    // add all neighbors indices around tile_index to tile_union
    pub fn _add_neighbors_to_tile_union(
        &self,
        tile_index: usize,
        tile_union: &mut [usize],
        n_near_tiles: &mut usize, // make sure we're passing in an indexable index for union (since capped at 27)
    ) {
        let tile = &self.tiles_struct._tiles[tile_index];
        for near_tile_idx in &tile.begin_tiles[0..tile.begin_len] {
            tile_union[*n_near_tiles] = *near_tile_idx;
            *n_near_tiles += 1;
        }
    }

    fn two_mut<T>(v: &mut [T], a: usize, b: usize) -> (&mut T, &mut T) {
        assert!(a != b);
        if a < b {
            let (left, right) = v.split_at_mut(b);
            (&mut left[a], &mut right[0])
        } else {
            let (left, right) = v.split_at_mut(a);
            (&mut right[0], &mut left[b])
        }
    }

    pub fn update_pair_nn(jet_a_idx: usize, jet_b_idx: usize, bj_jets: &mut [TiledJet]) {
        let (jet_a, jet_b) = ClusterSequence::two_mut(bj_jets, jet_a_idx, jet_b_idx);
        let dist = TiledJet::_bj_dist(jet_a, jet_b);

        let a_nn_dist = jet_a.nn_dist();
        let b_nn_dist = jet_b.nn_dist();

        if dist < a_nn_dist {
            jet_a.set_nn_dist(dist);
            jet_a.set_nn_jet(jet_b_idx);
        }

        if dist < b_nn_dist {
            jet_b.set_nn_dist(dist);
            jet_b.set_nn_jet(jet_a_idx);
        }
    }

    pub fn tiled_n2_cluster(&mut self) {
        self.initialize_tiles();
        let n = self.particles.len();
        let mut active: Vec<usize> = Vec::with_capacity(n);
        let mut bj_jets: Vec<TiledJet> = Vec::with_capacity(n); // stable arena
        let mut tile_union: Vec<usize> = vec![0; 3 * 9]; // giving vec length instead of just capacity

        // set bj info for all tiles
        for i in 0..n {
            let tile_jet_idx = self._tj_set_jetinfo(&mut bj_jets, None, i);
            bj_jets[tile_jet_idx].nn_jet_index = tile_jet_idx; // default nn is itself
            active.push(tile_jet_idx);
        }

        // set up all NN info
        self.tiles_struct._tiles.iter().for_each(|tile| {
            // set bj_dist
            let mut jet_a_next_idx = tile.head;

            while let Some(jet_a_idx) = jet_a_next_idx {
                let mut jet_b_next_idx = tile.head;
                while let Some(jet_b_idx) = jet_b_next_idx {
                    if jet_a_idx == jet_b_idx {
                        break;
                    }
                    ClusterSequence::update_pair_nn(jet_a_idx, jet_b_idx, &mut bj_jets);
                    jet_b_next_idx = bj_jets[jet_b_idx].next_jet;
                }
                jet_a_next_idx = bj_jets[jet_a_idx].next_jet;
            }
            // now we need to do it for the RH tiles since a NN might lie there for this tile
            // we do not check the LH tiles since these should already be computed by the prev
            // tile loop
            for rh_tile_idx in tile.rh_tiles.start..tile.rh_tiles.end {
                // rh tile should never be None
                let rh_tile = &self.tiles_struct._tiles[tile.begin_tiles[rh_tile_idx]];
                // do another loop over all tiles
                let mut jet_a_next_idx = tile.head;
                while let Some(jet_a_idx) = jet_a_next_idx {
                    // now we loop over all RH tile jets
                    let mut jet_b_next_idx = rh_tile.head;
                    while let Some(jet_b_idx) = jet_b_next_idx {
                        ClusterSequence::update_pair_nn(jet_a_idx, jet_b_idx, &mut bj_jets);
                        jet_b_next_idx = bj_jets[jet_b_idx].next_jet;
                    }
                    jet_a_next_idx = bj_jets[jet_a_idx].next_jet;
                }
            }
        });

        // now we want to create the DiJ distance metric where NN is J
        // this is easy since all the info is already computed prev
        // just loop over all jets
        // TODO fix this goofy ProxyJet shenanigans
        let mut di_j: Vec<f64> = vec![0.0; n];
        for pos in 0..active.len() {
            let id = active[pos];
            di_j[pos] = TiledJet::_bj_dij(&bj_jets[id], &bj_jets);
        }

        //start recombination loop
        for i in 0..n {
            let end_idx = n - 1 - i;
            let (min_idx, min_dij) = di_j[..=end_idx]
                .iter()
                .enumerate()
                .min_by_key(|(_, d)| d.to_bits())
                .unwrap();

            let jet_a_idx = active[min_idx];
            let jet_a_tile_idx = bj_jets[jet_a_idx].tile_index;
            let mut jet_b_idx = bj_jets[jet_a_idx].nn_jet_index;
            // normalize
            let min_dij = min_dij * self._invr2;

            // keep old pointer info after removing from tiles
            let mut old_b_tile_idx: Option<usize> = None;
            let mut old_b_jet_idx: Option<usize> = None;
            let mut new_b_jet_idx: Option<usize> = None;

            //check if jet_a has neighbor
            if jet_b_idx != jet_a_idx {
                let b_pos = active.iter().position(|&id| id == jet_b_idx).unwrap();

                // nn is the new index we want to store the recombined index at
                // modify jetB for BJ with new info
                let new_bj = self.do_jet_jet_recombination_step::<TiledJet>(
                    bj_jets[jet_a_idx].particle_index(),
                    bj_jets[jet_b_idx].particle_index(),
                    min_dij,
                    jet_b_idx,
                );

                // do jj now remove jetA from tiles
                // and remove old jetB
                // and add newly merged jetB
                self.bj_remove_from_tiles(&mut bj_jets, jet_a_idx);
                old_b_tile_idx = Some(bj_jets[jet_b_idx].tile_index);
                old_b_jet_idx = Some(jet_b_idx);
                self.bj_remove_from_tiles(&mut bj_jets, jet_b_idx);
                let part_idx = bj_jets[jet_b_idx].particle_index();

                let new_bj_rc_idx = self._tj_set_jetinfo(&mut bj_jets, Some(new_bj), part_idx);
                active[b_pos] = new_bj_rc_idx;
                new_b_jet_idx = Some(new_bj_rc_idx);
                jet_b_idx = new_bj_rc_idx;
            } else {
                self.do_jet_beam_recombination_step(bj_jets[jet_a_idx].particle_index(), min_dij);
                self.bj_remove_from_tiles(&mut bj_jets, jet_a_idx);
            };

            let mut n_near_tiles: usize = 0;
            // add neighbors to one search vector "tile_union"
            self._add_neighbors_to_tile_union(jet_a_tile_idx, &mut tile_union, &mut n_near_tiles);

            if let Some(new_b_idx) = new_b_jet_idx {
                // if we add other neighbors then add flag
                let mut need_sort = false;
                // if same tile then ignore
                if bj_jets[new_b_idx].tile_index != jet_a_tile_idx {
                    need_sort = true;
                    self._add_neighbors_to_tile_union(
                        bj_jets[new_b_idx].tile_index,
                        &mut tile_union,
                        &mut n_near_tiles,
                    );
                }

                //check if old jet needs to be added too
                if let Some(old_b_tile) = old_b_tile_idx
                    && old_b_tile != bj_jets[new_b_idx].tile_index
                    && old_b_tile != jet_a_tile_idx
                {
                    need_sort = true;
                    self._add_neighbors_to_tile_union(
                        old_b_tile,
                        &mut tile_union,
                        &mut n_near_tiles,
                    );
                }

                if need_sort {
                    let _ = &mut tile_union[0..n_near_tiles].sort();
                    //take sorted list and now compress repeated tiles
                    let mut compressed_idx = 1;
                    for idx in 1..n_near_tiles {
                        if tile_union[idx] != tile_union[compressed_idx - 1] {
                            tile_union[compressed_idx] = tile_union[idx];
                            compressed_idx += 1;
                        }
                    }
                    n_near_tiles = compressed_idx;
                }
            }

            //moves recombined min_idx with potential new NN jet at end_idx
            if min_idx != end_idx {
                active.swap(min_idx, end_idx);
                di_j.swap(min_idx, end_idx);
            }
            active.pop();
            di_j.pop();

            // now update new jetB NN by using tile union
            for &tile_union_neighbor in tile_union.iter().take(n_near_tiles) {
                // let tile = &self.tiles_struct._tiles[tile_union[tile_union_neighbor]];
                let tile = &self.tiles_struct._tiles[tile_union_neighbor];
                let mut jet_i_next_idx = tile.head;

                while let Some(jet_i_idx) = jet_i_next_idx {
                    // if jet_i had NN of jetA or jetB need to recalculate
                    if bj_jets[jet_i_idx].nn_jet_index == jet_a_idx
                        || old_b_jet_idx
                            .is_some_and(|old_b_id| bj_jets[jet_i_idx].nn_jet_index == old_b_id)
                    {
                        bj_jets[jet_i_idx].nn_dist = self.r2;
                        bj_jets[jet_i_idx].nn_jet_index = jet_i_idx;
                        // now find new NN in all tiles
                        for near_tile_idx in &tile.begin_tiles[0..tile.begin_len] {
                            // now go through bj over all jets in each tile
                            let mut jet_j_next_idx = self.tiles_struct._tiles[*near_tile_idx].head;
                            while let Some(jet_j_idx) = jet_j_next_idx {
                                let dist =
                                    TiledJet::_bj_dist(&bj_jets[jet_i_idx], &bj_jets[jet_j_idx]);
                                if dist < bj_jets[jet_i_idx].nn_dist() && jet_i_idx != jet_j_idx {
                                    bj_jets[jet_i_idx].nn_dist = dist;
                                    bj_jets[jet_i_idx].nn_jet_index = jet_j_idx;
                                }
                                jet_j_next_idx = bj_jets[jet_j_idx].next_jet;
                            }
                        }
                        if let Some(pos) = active.iter().position(|&id| id == jet_i_idx) {
                            di_j[pos] =
                                <TiledJet as ProxyJet>::_bj_dij(&bj_jets[jet_i_idx], &bj_jets);
                        }
                    }
                    // if JJ recomb check if new jet is closer than jeti's current NN
                    if let Some(new_b_idx) = new_b_jet_idx {
                        let dist = TiledJet::_bj_dist(&bj_jets[jet_i_idx], &bj_jets[jet_b_idx]);
                        if dist < bj_jets[jet_i_idx].nn_dist && jet_i_idx != new_b_idx {
                            bj_jets[jet_i_idx].nn_dist = dist;
                            bj_jets[jet_i_idx].nn_jet_index = jet_b_idx;
                            if let Some(pos) = active.iter().position(|&id| id == jet_i_idx) {
                                di_j[pos] =
                                    <TiledJet as ProxyJet>::_bj_dij(&bj_jets[jet_i_idx], &bj_jets);
                            }
                        }
                        if dist < bj_jets[new_b_idx].nn_dist && jet_i_idx != new_b_idx {
                            bj_jets[new_b_idx].nn_dist = dist;
                            bj_jets[new_b_idx].nn_jet_index = jet_i_idx;
                        }
                    }
                    jet_i_next_idx = bj_jets[jet_i_idx].next_jet;
                }
            }

            // update diJ for our new jetB
            // if jet_b_idx != min_idx {
            //     di_j[jet_b_idx] = TiledJet::_bj_dij(&bj_jets[jet_b_idx], &bj_jets);
            // }

            // // update all pointers of old tail tiles
            // let tiles = &self.tiles_struct._tiles[bj_jets[end_idx].tile_index];
            // for near_tiled_idx in &tiles.begin_tiles[0..tiles.begin_len] {
            //     let near_tile = &self.tiles_struct._tiles[*near_tiled_idx];
            //     let mut jet_j_next_idx = near_tile.head;
            //     while let Some(jet_j_idx) = jet_j_next_idx {
            //         if bj_jets[jet_j_idx].nn_jet_index == 0 {
            //             bj_jets[jet_j_idx].nn_jet_index = 0;
            //         }
            //         jet_j_next_idx = bj_jets[jet_j_idx].next_jet;
            //     }
            // }

            if let Some(new_b_idx) = new_b_jet_idx
                && let Some(pos) = active.iter().position(|&id| id == new_b_idx)
            {
                di_j[pos] = TiledJet::_bj_dij(&bj_jets[new_b_idx], &bj_jets);
            }
        }
    }

    #[inline(always)]
    #[cfg(feature = "simd")]
    fn set_nn_bj_vec_simd_and_dij(
        &self,
        rap_arr: &[f64],
        phi_arr: &[f64],
        nn_arr: &mut [u64],
        nn_dist_arr: &mut [f64],
        kt2_arr: &[f64],
        di_j: &mut [f64],
    ) {
        const LANES: usize = 16;
        let n = rap_arr.len();
        for i in 1..n {
            ClusterSequence::set_nn_crosscheck_simd::<LANES>(
                i,
                0..i - 1, // to s= i-1 inclusive
                rap_arr,
                phi_arr,
                self.r2,
                nn_dist_arr,
                nn_arr,
            );
        }

        ClusterSequence::_dij_simd::<LANES>(kt2_arr, nn_arr, nn_dist_arr, di_j);
    }

    #[inline(always)]
    #[cfg(feature = "simd")]
    pub fn set_nn_nocross_simd<const LANES: usize>(
        i: usize,
        range: std::ops::Range<usize>,
        rap_arr: &[f64],
        phi_arr: &[f64],
        r2: f64,
        nn_dist_arr: &mut [f64],
        nn_arr: &mut [u64],
    ) {
        let rap_i = rap_arr[i];
        let phi_i = phi_arr[i];
        let pi = Simd::<f64, LANES>::splat(std::f64::consts::PI);
        let two_pi = Simd::<f64, LANES>::splat(2.0 * std::f64::consts::PI);
        let inf = Simd::<f64, LANES>::splat(f64::INFINITY);
        let i_vec = Simd::<u64, LANES>::splat(i as u64);

        let mut best_d2 = Simd::<f64, LANES>::splat(r2);
        let mut best_j = Simd::<u64, LANES>::splat(i as u64);

        let lane_offsets = Simd::<u64, LANES>::from_array(std::array::from_fn(|k| k as u64));
        let mut j = range.start;

        while j + (LANES * 2) <= (range.end + 1) {
            let rap_j1 = Simd::<f64, LANES>::from_slice(&rap_arr[j..j + LANES]);
            let phi_j1 = Simd::<f64, LANES>::from_slice(&phi_arr[j..j + LANES]);
            let idx1 = Simd::<u64, LANES>::splat(j as u64) + lane_offsets;
            let dist1 =
                ClusterSequence::bj_dist_simd_direct(rap_i, phi_i, rap_j1, phi_j1, pi, two_pi);

            let mask1 = idx1.simd_ne(i_vec);
            let safe_dist1 = mask1.select(dist1, inf);

            let better1 = safe_dist1.simd_lt(best_d2);
            best_d2 = better1.select(safe_dist1, best_d2);
            best_j = better1.select(idx1, best_j);

            let j2 = j + LANES;
            let rap_j2 = Simd::<f64, LANES>::from_slice(&rap_arr[j2..j2 + LANES]);
            let phi_j2 = Simd::<f64, LANES>::from_slice(&phi_arr[j2..j2 + LANES]);
            let idx2 = Simd::<u64, LANES>::splat(j2 as u64) + lane_offsets;
            let dist2 =
                ClusterSequence::bj_dist_simd_direct(rap_i, phi_i, rap_j2, phi_j2, pi, two_pi);

            let mask2 = idx2.simd_ne(i_vec);
            let safe_dist2 = mask2.select(dist2, inf);

            let better2 = safe_dist2.simd_lt(best_d2);
            best_d2 = better2.select(safe_dist2, best_d2);
            best_j = better2.select(idx2, best_j);

            j += LANES * 2;
        }

        let mut nndist_min = best_d2.reduce_min();
        // Find the lane index for the winner
        let mask = best_d2.simd_eq(Simd::splat(nndist_min));
        let mut nn_min = mask.select(best_j, Simd::splat(u64::MAX)).reduce_min();

        // handle remaining elements (j..=end)
        while j <= range.end {
            if j != i {
                let drap = rap_i - rap_arr[j];
                let mut dphi = (phi_i - phi_arr[j]).abs();
                if dphi > std::f64::consts::PI {
                    dphi = 2.0 * std::f64::consts::PI - dphi;
                }
                let d2 = drap.mul_add(drap, dphi * dphi);

                if d2 < nndist_min {
                    nndist_min = d2;
                    nn_min = j as u64;
                }
            }
            j += 1;
        }

        nn_dist_arr[i] = nndist_min;
        nn_arr[i] = nn_min;
    }

    #[inline(always)]
    #[cfg(feature = "simd")]
    pub fn set_nn_crosscheck_simd<const LANES: usize>(
        i: usize,
        range: std::ops::Range<usize>,
        rap_arr: &[f64],
        phi_arr: &[f64],
        r2: f64,
        nndist_arr: &mut [f64],
        nn_arr: &mut [u64],
    ) {
        let rap_i = rap_arr[i];

        let phi_i = phi_arr[i];

        let mut j = range.start;

        let mut best_d2 = Simd::<f64, LANES>::splat(r2);
        let mut best_j = Simd::<u64, LANES>::splat(i as u64);

        // constant lane offsets [0,1,2,...]
        let lane_offsets: Simd<u64, LANES> = Simd::from_array(std::array::from_fn(|k| k as u64));

        let inf = Simd::<f64, LANES>::splat(f64::INFINITY);
        let i_u64 = i as u64;
        let i_vec = Simd::<u64, LANES>::splat(i_u64);

        let pi = Simd::<f64, LANES>::splat(std::f64::consts::PI);
        let two_pi = Simd::<f64, LANES>::splat(2.0 * std::f64::consts::PI);

        while j + LANES <= (range.end + 1) {
            let rap_j = Simd::<f64, LANES>::from_slice(&rap_arr[j..j + LANES]);
            let phi_j = Simd::<f64, LANES>::from_slice(&phi_arr[j..j + LANES]);
            let mut dist =
                ClusterSequence::bj_dist_simd_direct(rap_i, phi_i, rap_j, phi_j, pi, two_pi);

            let idx = Simd::<u64, LANES>::splat(j as u64) + lane_offsets;
            let is_self = idx.simd_eq(i_vec);
            dist = is_self.select(inf, dist);
            let better_for_i = dist.simd_lt(best_d2);
            best_d2 = better_for_i.select(dist, best_d2);
            best_j = better_for_i.select(idx, best_j);

            let old_d = Simd::<f64, LANES>::from_slice(&nndist_arr[j..j + LANES]);
            let improve = dist.simd_lt(old_d);

            let new_d = improve.select(dist, old_d);
            new_d.copy_to_slice(&mut nndist_arr[j..j + LANES]);

            // Update nn under same mask
            let old_nn = Simd::<u64, LANES>::from_slice(&nn_arr[j..j + LANES]);
            let i_vec = Simd::<u64, LANES>::splat(i as u64);
            let new_nn = improve.select(i_vec, old_nn);
            new_nn.copy_to_slice(&mut nn_arr[j..j + LANES]);

            j += LANES;
        }

        let mut nndist_min = r2;
        let mut nn_min = i as u64;

        while j <= range.end {
            if j == i {
                j += 1;
                continue;
            }
            let drap = rap_i - rap_arr[j];
            let mut dphi = (phi_i - phi_arr[j]).abs();
            if dphi > PI {
                dphi = 2.0 * PI - dphi
            }
            let d2 = drap.mul_add(drap, dphi * dphi);

            if d2 < nndist_min {
                nndist_min = d2;
                nn_min = j as u64;
            }

            if d2 < nndist_arr[j] {
                nndist_arr[j] = d2;
                nn_arr[j] = i as u64;
            }

            j += 1;
        }

        let bd = best_d2.to_array();
        let bj = best_j.to_array();
        for lane in 0..LANES {
            if bd[lane] < nndist_min {
                nndist_min = bd[lane];
                nn_min = bj[lane];
            }
        }

        nndist_arr[i] = nndist_min;
        nn_arr[i] = nn_min;
    }

    #[inline(always)]
    #[cfg(feature = "simd")]
    fn bj_dist_simd_direct<const LANES: usize>(
        rap_i: f64,
        phi_i: f64,
        rap_j: Simd<f64, LANES>,
        phi_j: Simd<f64, LANES>,
        pi_vec: Simd<f64, LANES>,
        two_pi_vec: Simd<f64, LANES>,
    ) -> Simd<f64, LANES> {
        let rap_i_vec = Simd::splat(rap_i);
        let phi_i_vec = Simd::splat(phi_i);

        let drap = rap_i_vec - rap_j;

        let mut dphi = (phi_i_vec - phi_j).abs();

        let mask = dphi.simd_gt(pi_vec);
        dphi = mask.select(two_pi_vec - dphi, dphi);

        dphi.mul_add(dphi, drap * drap)
    }

    #[inline(always)]
    #[cfg(feature = "simd")]
    fn _dij_simd<const LANES: usize>(
        kt2_arr: &[f64],
        nn_idx: &[u64],
        nn_dist: &[f64],
        di_j: &mut [f64],
    ) {
        let n = kt2_arr.len();
        let mut j = 0;

        while j + LANES <= n {
            let kt2_i = Simd::<f64, LANES>::from_slice(&kt2_arr[j..j + LANES]);
            let dist_i = Simd::<f64, LANES>::from_slice(&nn_dist[j..j + LANES]);
            let indices = Simd::<u64, LANES>::from_slice(&nn_idx[j..j + LANES]);

            let indices_usize = indices.cast::<usize>();

            let kt2_b = Simd::<f64, LANES>::gather_or(kt2_arr, indices_usize, kt2_i);

            let kt2_min = kt2_i.simd_min(kt2_b);
            let res = kt2_min * dist_i;

            res.copy_to_slice(&mut di_j[j..j + LANES]);

            j += LANES;
        }

        for i in j..n {
            let idx = nn_idx[i] as usize;
            di_j[i] = nn_dist[i] * kt2_arr[i].min(kt2_arr[idx]);
        }
    }

    #[inline(always)]
    fn set_nn_and_dij<J: ProxyJet>(&mut self, bj_jets: &mut [J], di_j: &mut [f64]) {
        let n = bj_jets.len();
        for i in 1..n {
            self.bj_set_nn_crosscheck(&mut bj_jets[0..=i]);
        }

        di_j.iter_mut().enumerate().for_each(|(i, jet)| {
            *jet = J::_bj_dij(&bj_jets[i], bj_jets);
        });
    }

    #[inline(always)]
    #[cfg(feature = "simd")]
    fn update_neighbors_recombined_full_simd<const LANES: usize>(
        &self,
        ctx: &mut ClusterContext,
        m_idx: usize,
        swapped_idx: usize,
        old_tail_idx: usize,
        n_active: usize,
    ) {
        let ClusterContext {
            rap,
            phi,
            kt2,
            nn_dist,
            di_j,
            nn_idx,
        } = ctx;

        let rap_m_val = rap[m_idx];
        let phi_m_val = phi[m_idx];
        let rap_m = Simd::<f64, LANES>::splat(rap_m_val);
        let phi_m = Simd::<f64, LANES>::splat(phi_m_val);
        let pi = Simd::<f64, LANES>::splat(std::f64::consts::PI);
        let two_pi = Simd::<f64, LANES>::splat(2.0 * std::f64::consts::PI);
        let lane_offsets = Simd::<u64, LANES>::from_array(std::array::from_fn(|i| i as u64));

        let mut m_best_d2 = Simd::<f64, LANES>::splat(self.r2);
        let mut m_best_nn = Simd::<u64, LANES>::splat(m_idx as u64);

        // the jets that lost a neighbor with the merge shouldnt be too high
        let mut needs_rescan = Vec::with_capacity(8);

        let mut j = 0;
        while j + LANES <= n_active {
            let rap_j = Simd::<f64, LANES>::from_slice(&rap[j..j + LANES]);
            let phi_j = Simd::<f64, LANES>::from_slice(&phi[j..j + LANES]);
            let idx_j = Simd::<u64, LANES>::splat(j as u64) + lane_offsets;

            let drap = rap_m - rap_j;
            let mut dphi = (phi_m - phi_j).abs();
            dphi = dphi.simd_gt(pi).select(two_pi - dphi, dphi);
            let d2 = dphi.mul_add(dphi, drap * drap);

            let mask_not_m = idx_j.simd_ne(Simd::<u64, LANES>::splat(m_idx as u64));
            let safe_d2 = mask_not_m.select(d2, Simd::splat(f64::INFINITY));

            let better_for_m = safe_d2.simd_lt(m_best_d2);
            m_best_d2 = better_for_m.select(safe_d2, m_best_d2);
            m_best_nn = better_for_m.select(idx_j, m_best_nn);

            let d2_vals = safe_d2.to_array();
            for (lane, d2) in d2_vals.iter().enumerate().take(LANES) {
                let i = j + lane;
                if i == m_idx || i >= n_active {
                    continue;
                }

                let cur_nn = nn_idx[i] as usize;
                if cur_nn == swapped_idx || cur_nn == old_tail_idx {
                    needs_rescan.push((i, *d2));
                } else if *d2 < nn_dist[i] {
                    // if m is the lowest distance then no need to find N as other have not changed
                    nn_dist[i] = *d2;
                    nn_idx[i] = m_idx as u64;
                    di_j[i] = nn_dist[i] * kt2[i].min(kt2[m_idx]);
                }
            }
            j += LANES;
        }

        let mut final_m_d2 = m_best_d2.reduce_min();
        let mask = m_best_d2.simd_eq(Simd::splat(final_m_d2));
        let mut final_m_nn = mask.select(m_best_nn, Simd::splat(u64::MAX)).reduce_min();

        while j < n_active {
            if j != m_idx {
                let drap = rap_m_val - rap[j];
                let mut dphi = (phi_m_val - phi[j]).abs();
                if dphi > std::f64::consts::PI {
                    dphi = 2.0 * std::f64::consts::PI - dphi;
                }
                let d2 = drap.mul_add(drap, dphi * dphi);

                if d2 < final_m_d2 {
                    final_m_d2 = d2;
                    final_m_nn = j as u64;
                }

                let cur_nn = nn_idx[j] as usize;
                if cur_nn == swapped_idx || cur_nn == old_tail_idx {
                    needs_rescan.push((j, d2));
                } else if d2 < nn_dist[j] {
                    nn_dist[j] = d2;
                    nn_idx[j] = m_idx as u64;
                    di_j[j] = d2 * kt2[j].min(kt2[m_idx]);
                }
            }
            j += 1;
        }

        nn_dist[m_idx] = final_m_d2;
        nn_idx[m_idx] = final_m_nn;
        di_j[m_idx] = final_m_d2 * kt2[m_idx].min(kt2[final_m_nn as usize]);

        for (idx, dist_to_m) in needs_rescan {
            ClusterSequence::set_nn_nocross_simd::<LANES>(
                idx,
                0..n_active - 1,
                rap,
                phi,
                self.r2,
                nn_dist,
                nn_idx,
            );

            if dist_to_m < nn_dist[idx] {
                nn_dist[idx] = dist_to_m;
                nn_idx[idx] = m_idx as u64;
            }

            let final_nn = nn_idx[idx] as usize;
            di_j[idx] = nn_dist[idx] * kt2[idx].min(kt2[final_nn]);
        }
    }

    #[inline(always)]
    #[cfg(feature = "simd")]
    fn update_neighbors_beam_simd<const LANES: usize>(
        &self,
        ctx: &mut ClusterContext,
        swapped_idx: usize,
        old_tail_idx: usize,
        n_active: usize,
    ) {
        let ClusterContext {
            rap,
            phi,
            kt2,
            nn_dist,
            di_j,
            nn_idx,
        } = ctx;

        for i in 0..n_active {
            let needs_recompute =
                nn_idx[i] == swapped_idx as u64 || nn_idx[i] == old_tail_idx as u64;

            if needs_recompute {
                ClusterSequence::set_nn_nocross_simd::<LANES>(
                    i,
                    0..n_active - 1,
                    rap,
                    phi,
                    self.r2,
                    nn_dist,
                    nn_idx,
                );

                // Update di_j for this jet now that it has a new neighbor
                let new_nn = nn_idx[i] as usize;
                di_j[i] = nn_dist[i] * kt2[i].min(kt2[new_nn]);
            }
        }
    }

    #[cfg(feature = "simd")]
    pub fn simple_n2_cluster_simd(&mut self) {
        //init all arrays instead of setting bj_set_jet info

        const LANES: usize = 16;
        let n = self.particles.len();
        let mut particles_index: Vec<usize> = (0..n).collect();
        let mut buffer = ClusterBuffer::new(n);
        let mut ctx = buffer.get_context();

        let r2 = self.r2;
        let algo = &self.jetdef.algorithm;
        let extra = self.jetdef.extra_param;

        ctx.nn_dist.fill(r2);

        assert_eq!(ctx.rap.len(), n);
        assert_eq!(ctx.phi.len(), n);
        assert_eq!(ctx.nn_idx.len(), n);
        assert_eq!(ctx.kt2.len(), n);

        for i in 0..n {
            let jet = &self.particles[i];

            ctx.rap[i] = jet.rap();
            ctx.phi[i] = jet.phi();
            ctx.nn_idx[i] = i as u64;

            ctx.kt2[i] = ClusterSequence::jet_scale_for_algorithm(jet, algo, extra);
        }

        self.set_nn_bj_vec_simd_and_dij(
            ctx.rap,
            ctx.phi,
            ctx.nn_idx,
            ctx.nn_dist,
            ctx.kt2,
            ctx.di_j,
        );

        // begin looping thru entire vector

        for i in 0..n {
            //find min distance
            let end_idx = n - 1 - i;
            let (mut min_idx, min_dij) = ctx
                .di_j
                .iter()
                .take(end_idx + 1)
                .enumerate()
                .min_by_key(|a| a.1.to_bits())
                .unwrap();

            let mut jet_b_idx = ctx.nn_idx[min_idx] as usize;

            // normalize
            let min_dij = min_dij * self._invr2;

            //check if jet_a has neighbor
            if jet_b_idx != min_idx {
                if min_idx < jet_b_idx {
                    std::mem::swap(&mut min_idx, &mut jet_b_idx);
                }
                // nn is the new index we want to store the recombined index at
                // modify jetB for BJ with new info

                let jet_i = &self.particles[particles_index[min_idx]];
                let jet_j = &self.particles[particles_index[jet_b_idx]];

                let new_idx = self.particles.len();
                let hist_a_idx = jet_i.cluster_hist_index();
                let hist_b_idx = jet_j.cluster_hist_index();

                let new_jet = self.jetdef.recombine(jet_i, jet_j);

                //update our arrays for new_jet
                ctx.rap[jet_b_idx] = new_jet.rap();
                ctx.phi[jet_b_idx] = new_jet.phi();
                ctx.nn_idx[jet_b_idx] = jet_b_idx as u64; //all nn will be set to None aka itself
                ctx.nn_dist[jet_b_idx] = self.r2;
                ctx.kt2[jet_b_idx] = ClusterSequence::jet_scale_for_algorithm(
                    &new_jet,
                    &self.jetdef.algorithm,
                    self.jetdef.extra_param,
                );
                particles_index[jet_b_idx] = new_idx;

                self.particles.push(new_jet);
                self.particles[new_idx].set_cluster_hist_index(self._history.len());

                self.add_recomb_to_history(
                    std::cmp::min(hist_a_idx, hist_b_idx),
                    std::cmp::max(hist_a_idx, hist_b_idx),
                    new_idx,
                    min_dij,
                );
            } else {
                self.add_recomb_to_history(
                    self.particles[particles_index[min_idx]].cluster_hist_index(),
                    BEAMJET,
                    INVALID,
                    min_dij,
                );
            };

            if min_idx != end_idx {
                ctx.rap[min_idx] = ctx.rap[end_idx];
                ctx.phi[min_idx] = ctx.phi[end_idx];
                ctx.kt2[min_idx] = ctx.kt2[end_idx];
                ctx.nn_dist[min_idx] = ctx.nn_dist[end_idx];
                ctx.nn_idx[min_idx] = ctx.nn_idx[end_idx];
                ctx.di_j[min_idx] = ctx.di_j[end_idx];
                particles_index[min_idx] = particles_index[end_idx];
            }

            let tail = end_idx;

            if jet_b_idx != min_idx {
                self.update_neighbors_recombined_full_simd::<LANES>(
                    &mut ctx, jet_b_idx, min_idx, tail, tail,
                );
            } else {
                self.update_neighbors_beam_simd::<LANES>(&mut ctx, min_idx, tail, tail);
            }
        }
    }

    pub fn simple_n2_cluster<J: ProxyJet>(&mut self) {
        let n = self.particles.len();

        let mut bj_jets: Vec<J> = Vec::with_capacity(n);

        // set brief jet info
        // TODO; prolly will change when multithreading because of iterator and vec push
        self.particles.iter().enumerate().for_each(|(i, jet)| {
            // the idx of NN will now be itself
            let bj: J = ProxyJet::_bj_set_jetinfo(
                i,
                jet,
                self.r2,
                ClusterSequence::jet_scale_for_algorithm(
                    jet,
                    &self.jetdef.algorithm,
                    self.jetdef.extra_param,
                ),
                i,
            );
            bj_jets.push(bj);
        });

        let mut di_j: Vec<f64> = vec![0.0; n];
        self.set_nn_and_dij(&mut bj_jets, &mut di_j);

        for i in 0..n {
            //find min distance
            let end_idx = n - 1 - i;
            let (mut min_idx, min_dij) = di_j
                .iter()
                .take(end_idx + 1)
                .enumerate()
                .min_by_key(|a| a.1.to_bits())
                .unwrap();

            // assert_eq!(min_i, min_idx);
            // assert_eq!(min_dij, *min_di);

            // let jet_a = &bj_jets[min_idx];
            // let mut jet_b_idx = jet_a.nn_jet_index(); //THIS IS AN OPTION

            // normalize
            let min_dij = min_dij * self._invr2;
            let jet_a = &bj_jets[min_idx];
            let mut jet_b_idx = jet_a.nn_jet_index();

            //check if jet_a has neighbor
            if jet_b_idx != min_idx {
                if min_idx < jet_b_idx {
                    std::mem::swap(&mut min_idx, &mut jet_b_idx);
                }
                // nn is the new index we want to store the recombined index at
                // modify jetB for BJ with new info

                let new_bj = self.do_jet_jet_recombination_step(
                    bj_jets[min_idx].particle_index(),
                    bj_jets[jet_b_idx].particle_index(),
                    min_dij,
                    jet_b_idx,
                );
                bj_jets[jet_b_idx] = new_bj;
            } else {
                self.do_jet_beam_recombination_step(bj_jets[min_idx].particle_index(), min_dij);
            };

            //moves recombined min_idx with potential new NN jet at end_idx
            //copy tail values to jetA
            bj_jets.swap(min_idx, end_idx);
            di_j.swap(min_idx, end_idx);

            //check entire list and update NN of references to new jet
            //very slow O(n) each time, compared to implicit pointer copying in C++
            // for jet in bj_jets.iter_mut().take(end_idx) {
            //     if jet.nn_jet_index() == end_idx {
            //         jet.set_nn_jet(min_idx);
            //     }
            // }

            // find all NN of recombined jet_b
            let tail = end_idx;
            if jet_b_idx != min_idx {
                // use these vecs to prove assertion bounds
                let jets = &mut bj_jets[..tail];
                let dijs = &mut di_j[..tail];

                assert!(jet_b_idx < jets.len());

                for i in 0..jet_b_idx {
                    let nn_idx = jets[i].nn_jet_index();

                    if nn_idx == min_idx || nn_idx == jet_b_idx {
                        let (new_nn, new_dist) = self.bj_set_nn_nocross(i, 0, tail, jets);

                        let jet_i = &mut jets[i];
                        jet_i.set_nn_jet(new_nn);
                        jet_i.set_nn_dist(new_dist);

                        dijs[i] = J::_bj_dij(&jets[i], jets);
                    }

                    let jet_ib_dist = J::_bj_dist(&jets[i], &jets[jet_b_idx]);

                    let jet_i = &mut jets[i];

                    if jet_i.nn_jet_index() == tail {
                        jet_i.set_nn_jet(min_idx);
                    }

                    if jet_ib_dist < jet_i.nn_dist() {
                        jet_i.set_nn_dist(jet_ib_dist);
                        jet_i.set_nn_jet(jet_b_idx);
                        dijs[i] = J::_bj_dij(&jets[i], jets);
                    }

                    let jet_b = &mut jets[jet_b_idx];
                    if jet_ib_dist < jet_b.nn_dist() {
                        jet_b.set_nn_dist(jet_ib_dist);
                        jet_b.set_nn_jet(i);
                    }
                }

                if jets[jet_b_idx].nn_jet_index() == tail {
                    jets[jet_b_idx].set_nn_jet(min_idx);
                }

                // LLVM also proves `i < tail` (since the loop goes to `tail`),
                // and `tail == jets.len()`. Again, zero bounds checks inside the loop.
                for i in jet_b_idx + 1..tail {
                    let nn_idx = jets[i].nn_jet_index();

                    if nn_idx == min_idx || nn_idx == jet_b_idx {
                        let (new_nn, new_dist) = self.bj_set_nn_nocross(i, 0, tail, jets);

                        let jet_i = &mut jets[i];
                        jet_i.set_nn_jet(new_nn);
                        jet_i.set_nn_dist(new_dist);

                        dijs[i] = J::_bj_dij(&jets[i], jets);
                    }

                    let jet_ib_dist = J::_bj_dist(&jets[i], &jets[jet_b_idx]);

                    let jet_i = &mut jets[i];

                    if jet_i.nn_jet_index() == tail {
                        jet_i.set_nn_jet(min_idx);
                    }

                    if jet_ib_dist < jet_i.nn_dist() {
                        jet_i.set_nn_dist(jet_ib_dist);
                        jet_i.set_nn_jet(jet_b_idx);
                        dijs[i] = J::_bj_dij(&jets[i], jets);
                    }

                    let jet_b = &mut jets[jet_b_idx];
                    if jet_ib_dist < jet_b.nn_dist() {
                        jet_b.set_nn_dist(jet_ib_dist);
                        jet_b.set_nn_jet(i);
                    }
                }

                dijs[jet_b_idx] = J::_bj_dij(&jets[jet_b_idx], jets);
            } else {
                //check if old JetA was NN to jet being mapped
                let jets = &mut bj_jets[..tail];
                let dijs = &mut di_j[..tail];

                for i in 0..tail {
                    let jet_i_nn_idx = jets[i].nn_jet_index();

                    if jet_i_nn_idx == min_idx {
                        let (new_nn, new_dist) = self.bj_set_nn_nocross(i, 0, tail, jets);

                        let jet_i = &mut jets[i];
                        jet_i.set_nn_jet(new_nn);
                        jet_i.set_nn_dist(new_dist);

                        dijs[i] = J::_bj_dij(&jets[i], jets);
                    } else if jet_i_nn_idx == tail {
                        jets[i].set_nn_jet(min_idx);
                    }
                }
            };
        }
    }

    fn do_jet_jet_recombination_step<J: ProxyJet>(
        &mut self,
        jet_a_idx: usize,
        jet_b_idx: usize,
        min_dij: f64,
        bj_jet_b_idx: usize,
    ) -> J {
        let new_jet = self
            .jetdef
            .recombine(&self.particles[jet_a_idx], &self.particles[jet_b_idx]);

        let new_bj: J = ProxyJet::_bj_set_jetinfo(
            self.particles.len(),
            &new_jet,
            self.r2,
            ClusterSequence::jet_scale_for_algorithm(
                &new_jet,
                &self.jetdef.algorithm,
                self.jetdef.extra_param,
            ),
            bj_jet_b_idx,
        );

        self.particles.push(new_jet);

        let new_idx = self.particles.len() - 1;
        let hist_a_idx = self.particles[jet_a_idx].cluster_hist_index();
        let hist_b_idx = self.particles[jet_b_idx].cluster_hist_index();

        self.particles[new_idx].set_cluster_hist_index(self._history.len());

        self.add_recomb_to_history(
            std::cmp::min(hist_a_idx, hist_b_idx),
            std::cmp::max(hist_a_idx, hist_b_idx),
            new_idx,
            min_dij,
        );

        new_bj
    }

    fn do_jet_beam_recombination_step(&mut self, jet_index: usize, d_ij: f64) {
        self.add_recomb_to_history(
            self.particles[jet_index].cluster_hist_index(),
            BEAMJET,
            INVALID,
            d_ij,
        );
    }

    fn add_recomb_to_history(
        &mut self,
        parent1: usize,
        parent2: usize,
        jet_index: usize,
        d_ij: f64,
    ) {
        let element: HistoryElement = HistoryElement {
            parent1,
            parent2,
            jet_index,
            d_ij,
            child: INVALID,
            curr_max_d_ij: f64::max(d_ij, self._history.last().unwrap().curr_max_d_ij), //TODO: verify options and stuff
        };

        self._history.push(element);

        let child_hist_index: usize = self._history.len() - 1;
        if self._history[parent1].child != INVALID {
            panic!(
                "Child {} history index for {} already set for parent 1",
                self._history[parent1].child, parent1
            );
        }
        self._history[parent1].child = child_hist_index;

        if !(parent2 == BEAMJET || parent2 == INVALID || parent2 == INEXISTENT_PARENT) {
            if self._history[parent2].child != INVALID {
                panic!(
                    "Child {} history index for {} already set for parent 2",
                    self._history[parent2].child, parent2
                );
            }
            self._history[parent2].child = child_hist_index;
        }

        if jet_index != INVALID {
            if jet_index == BEAMJET || jet_index == INEXISTENT_PARENT {
                panic!("Jet_index has incorrect index (Beam jet or Inexistent parent)");
            }
            self.particles[jet_index].set_cluster_hist_index(child_hist_index);
        }
    }

    // -----------STATIC METHODS USING PROXYJET--------------

    #[inline]
    fn bj_set_nn_nocross<J: ProxyJet>(
        &self,
        curr_idx: usize,
        head_idx: usize,
        tail_idx: usize,
        jets: &[J],
    ) -> (usize, f64) {
        let mut nn_dist = self.r2;
        let mut nn: usize = curr_idx;

        // attempt to remove bound checks
        let len = jets.len();
        assert!(curr_idx < len);
        assert!(tail_idx <= len);

        // this part is immutable only change at the end
        let jets_ref = &jets;
        let curr_jet = &jets_ref[curr_idx];

        for jet_b_idx in head_idx..curr_idx {
            let dist = J::_bj_dist(curr_jet, &jets_ref[jet_b_idx]);
            if dist < nn_dist {
                nn_dist = dist;
                nn = jet_b_idx;
            }
        }

        for jet_b_idx in curr_idx + 1..tail_idx {
            let dist = J::_bj_dist(curr_jet, &jets_ref[jet_b_idx]);
            if dist < nn_dist {
                nn_dist = dist;
                nn = jet_b_idx;
            }
        }
        (nn, nn_dist)
    }

    #[inline]
    fn bj_set_nn_crosscheck<J: ProxyJet>(&self, jets: &mut [J]) {
        let mut nn_dist = self.r2;
        let n = jets.len() - 1;
        let mut nn: usize = n;
        for i in 0..n {
            let dist = J::_bj_dist(&jets[i], &jets[n]);
            if dist < nn_dist {
                nn_dist = dist;
                nn = i;
            }
            if dist < jets[i].nn_dist() {
                jets[i].set_nn_dist(dist);
                jets[i].set_nn_jet(n);
            }
        }
        jets[n].set_nn_dist(nn_dist);
        jets[n].set_nn_jet(nn);
    }
}
