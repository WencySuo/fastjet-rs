use log::warn;
use std::cell::RefCell;
use std::panic;
use std::rc::Rc;

use crate::constants::PI;
use crate::constants::TWO_PI;
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

// TODO: implement other algorithms
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
        let mut kt2 = *jet.kt2();
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
                //TODO: add some logic here to determine which jettype to use
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
                        if *jet.perp2() >= dcut {
                            Some(jet)
                        } else {
                            None
                        }
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
                    if *jet.perp2() >= dcut {
                        Some(jet)
                    } else {
                        None
                    }
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
            let rap = *particle.rap();
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
        println!("busiest_bin: {}", busiest_bin);

        //more magic numbers
        let allowed_max_fraction: f64 = 0.25;
        let min_multiplicity: f64 = 4.0;

        //find how much we can stuff in edge bins
        let mut allowed_max_cumul = min_multiplicity
            .max(*busiest_bin as f64 * allowed_max_fraction)
            .floor();

        println!("allowed_max_cumul: {}", allowed_max_cumul);

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
        println!("after left scan cumul_2: {}", cumul_2);

        let ibin_lo = ubin;

        let mut cum_hi = 0.0;
        // do same as right except from hi to lo
        ubin = 0;
        for i in (0..nbins).rev() {
            cum_hi += counts[i] as f64;
            if cum_hi >= allowed_max_cumul {
                //RHS side of last bin
                let y = (i as f64 - nrap as f64 + 1.0) as f64;
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
        println!(
            "n_tiles_phi: {}, _tile_size_phi: {}",
            self.tiles_struct._n_tiles_phi, self.tiles_struct._tile_size_phi
        );

        // find the min and max rap/eta that we should be using for this analysis
        let (min_rap, max_rap, _cumul_2) = self.min_max_rap();
        println!(
            "min_rap: {}, max_rap: {}, _cumul_2: {}",
            min_rap, max_rap, _cumul_2
        );

        //find min/max values for these tiles and figure out why _tiles_eta_min isnt just min_rap?
        self.tiles_struct._tiles_ieta_min =
            (min_rap / self.tiles_struct._tile_size_eta).floor() as isize;
        self.tiles_struct._tiles_ieta_max =
            (max_rap / self.tiles_struct._tile_size_eta).floor() as isize;
        self.tiles_struct._tiles_eta_min =
            self.tiles_struct._tiles_ieta_min as f64 * self.tiles_struct._tile_size_eta;
        self.tiles_struct._tiles_eta_max =
            self.tiles_struct._tiles_ieta_max as f64 * self.tiles_struct._tile_size_eta;

        println!("_tiles_ieta_min={}", self.tiles_struct._tiles_ieta_min);
        println!("_tiles_ieta_max={}", self.tiles_struct._tiles_ieta_max);
        println!("_tiles_eta_min={}", self.tiles_struct._tiles_eta_min);
        println!("_tiles_eta_max={}", self.tiles_struct._tiles_eta_max);

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
        for ieta in self.tiles_struct._tiles_ieta_min..self.tiles_struct._tiles_ieta_max {
            for iphi in 0..self.tiles_struct._n_tiles_phi {
                let mut tile_idx = 0;
                let tile: &mut Tile = &mut tiles[self._tile_int_index(ieta, iphi)];

                let print_tile_idx = self._tile_int_index(ieta, iphi);
                println!(
                    "starting cross-referencing at ieta {} and iphi {} for tile {}",
                    ieta, iphi, print_tile_idx
                );

                // first tile has no HEAD
                tile.head = Option::None;

                //begin_tiles ppoints to neighboring tiles including itself?
                tile.begin_tiles[tile_idx] = self._tile_int_index(ieta, iphi);

                // the surrounding tiles excludes self
                // can just use slices in rust
                tile_idx += 1;
                //tile.surrounding_tiles = &tile.begin_tiles[tile_idx..];

                // after first loop of ieta we can now check prev column
                // now left right middle
                if ieta > self.tiles_struct._tiles_ieta_min {
                    for idphi in [-1isize, 0, 1] {
                        tile_idx += 1;
                        tile.begin_tiles[tile_idx] = self._tile_int_index(ieta - 1, iphi + idphi);
                    }
                }

                // note: phi doesnt need bound checks since it mod 2pi
                // LHS includes under curr elem
                tile_idx += 1;
                tile.begin_tiles[tile_idx] = self._tile_int_index(ieta, iphi - 1);

                // first RHS is above curr elem
                tile_idx += 1;
                tile.begin_tiles[tile_idx] = self._tile_int_index(ieta, iphi + 1);
                println!(
                    "setup first R index={}",
                    self._tile_int_index(ieta, iphi + 1)
                );

                // only set last R if we are not at max
                if ieta < self.tiles_struct._tiles_ieta_max {
                    for idphi in [-1isize, 0, 1] {
                        tile_idx += 1;
                        tile_idx = tile_idx % 9 as usize;
                        let print_tile_idx = self._tile_int_index(ieta + 1, iphi + idphi);
                        println!("print out of bounds tile index: {}", print_tile_idx);
                        tile.begin_tiles[tile_idx] = self._tile_int_index(ieta + 1, iphi + idphi);
                        // println!(
                        //     "setup remaining R's if ieta is not max {}",
                        //     self._tile_int_index(ieta + 1, iphi + idphi)
                        // );
                    }
                }

                // TODO: since i think this is unitialized in prev config , aka at correct points
                // check if this does not have runtime errors?
                // TODO: this is not fully defined to be these constants since need to accound for when
                // ieta > self._tiles_ieta_min  and when ieta < self._tiles_ieta_max
                tile.surrounding_tiles = 1..9;
                tile.rh_tiles = 5..9;
            }
        }
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
        tiled_jet: Option<TiledJet>,
        particle_index: usize,
    ) -> Rc<RefCell<TiledJet>> {
        let jet = if let Some(tj) = tiled_jet {
            Rc::new(RefCell::new(tj))
        } else {
            let pseudo_jet = &self.particles[particle_index];
            let new_tj = TiledJet::_bj_set_jetinfo(
                particle_index,
                pseudo_jet,
                self.r2,
                ClusterSequence::jet_scale_for_algorithm(
                    pseudo_jet,
                    &self.jetdef.algorithm,
                    self.jetdef.extra_param,
                ),
            );
            Rc::new(RefCell::new(new_tj))
        };

        let (eta, phi) = {
            let j = jet.borrow();
            (j.eta(), j.phi())
        };

        jet.borrow_mut().tile_index = self._tile_index(eta, phi);
        let tile_index = jet.borrow().tile_index;

        // need to add jet pointer to linked List of jets
        {
            let mut j = jet.borrow_mut();
            j.prev_jet = None;
            j.next_jet = self.tiles_struct._tiles[tile_index].head.clone();
        }

        if let Some(next) = self.tiles_struct._tiles[tile_index].head.as_ref() {
            next.borrow_mut().prev_jet = Some(Rc::clone(&jet));
        }

        self.tiles_struct._tiles[tile_index].head = Some(Rc::clone(&jet));
        jet
    }

    pub fn bj_remove_from_tiles(&mut self, jet: &Rc<RefCell<TiledJet>>) {
        let tile = &mut self.tiles_struct._tiles[jet.borrow().tile_index];

        // if prev is null then jet is head
        // TODO: check if we really need to clone in both if lets
        if let Some(ref prev) = jet.borrow().prev_jet {
            // link prev pointer to next pointer to remove jet
            prev.borrow_mut().next_jet = jet.borrow().next_jet.clone();
        } else {
            tile.head = jet.borrow().next_jet.clone();
        }

        if let Some(ref next) = jet.borrow().next_jet {
            next.borrow_mut().prev_jet = jet.borrow().prev_jet.clone();
        }

        //TODO check if we need to drop jet that we just added?
    }

    // add all neighbors indices around tile_index to tile_union
    pub fn _add_neighbors_to_tile_union(
        &mut self,
        tile_index: usize,
        tile_union: &mut [usize],
        n_near_tiles: &mut usize,
    ) {
        for near_tile_idx in self.tiles_struct._tiles[tile_index].begin_tiles {
            tile_union[*n_near_tiles] = near_tile_idx;
            *n_near_tiles += 1;
        }
    }

    pub fn tiled_n2_cluster(&mut self) {
        //now we follow similar patterns to simple_n2_cluster
        self.initialize_tiles();

        let n = self.particles.len();

        let mut bj_jets: Vec<Rc<RefCell<TiledJet>>> = Vec::with_capacity(n);

        // TODO find out what tile union does lmao
        // TODO set constant for n_tile neighbors which is NINE
        let mut tile_union: Vec<usize> = Vec::with_capacity(3 * 9);

        // Move tiles out so we don't hold a &mut borrow into self while calling methods on self
        // TODO: there is some funky self borrow checker stuff when setting info
        // because moving entire struct this should be solved once we move things like tile_index
        // into its own class relying only on ClusterSequenceTiles
        let tiles = std::mem::take(&mut self.tiles_struct._tiles);
        // let particles = &self.particles;

        //set bj info for all tiles
        for i in 0..n {
            let tile_jet = self._tj_set_jetinfo(None, i);
            tile_jet.borrow_mut().bj_jet_index = i;
            bj_jets.push(tile_jet);
        }
        // for (i, pseudo_jet) in particles.iter().enumerate() {
        //     let mut tile_jet = self._tj_set_jetinfo(None, i);
        //     // reluctantly add jet index for getting index in linked lists
        //     tile_jet.borrow_mut().bj_jet_index = i;
        //     bj_jets.push(tile_jet);
        // }

        self.tiles_struct._tiles = tiles;

        // now set up all NN info, looping through tiles and their neighbors

        self.tiles_struct._tiles.iter().for_each(|tile| {
            // N2 loop over all NN setting their bj_dist
            // cloning here is proper since its just an Rc but most likely hurt perf
            let mut jet_a_next = tile.head.clone();
            while let Some(jet_a) = jet_a_next {
                //jet_a cannot be null
                let mut jet_b_next = tile.head.clone();
                while let Some(jet_b) = jet_b_next
                    && Rc::ptr_eq(&jet_a, &jet_b)
                {
                    let dist = TiledJet::_bj_dist(&jet_a.borrow(), &jet_b.borrow());
                    if dist < jet_a.borrow().nn_dist {
                        jet_a.borrow_mut().nn_dist = dist;
                        jet_a.borrow_mut().nn_jet_index = jet_b.borrow().nn_jet_index;
                    }
                    if dist < jet_b.borrow().nn_dist {
                        jet_b.borrow_mut().nn_dist = dist;
                        jet_b.borrow_mut().nn_jet_index = jet_a.borrow().nn_jet_index;
                    }
                    jet_b_next = jet_a.borrow().next_jet.clone();
                }
                // cloning here is proper since its just an Rc but also slow
                jet_a_next = jet_a.borrow().next_jet.clone();
            }
            // now we need to do it for the RH tiles since a NN might lie there for this tile
            // we do not check the LH tiles since these should already be computed by the prev
            // tile loop
            //TODO rh tiles may not be range in future impl on beginning tiles
            for rh_tile_idx in tile.rh_tiles.start..tile.rh_tiles.end {
                // rh tile should never be None
                let rh_tile = &self.tiles_struct._tiles[rh_tile_idx];
                // do another loop over all tiles
                let mut jet_a_next = tile.head.clone();
                while let Some(jet_a) = jet_a_next {
                    // now we loop over all RH tile jets
                    let mut jet_b_next = rh_tile.head.clone();
                    while let Some(jet_b) = jet_b_next {
                        // TODO add _bj_dist to tiledJet impl placeholder for now
                        // BriefJet::_bj_dist(jet_a.borrow(), jet_b.borrow());
                        let dist = TiledJet::_bj_dist(&jet_a.borrow(), &jet_b.borrow());
                        if dist < jet_a.borrow().nn_dist {
                            jet_a.borrow_mut().nn_dist = dist;
                            jet_a.borrow_mut().nn_jet_index = jet_b.borrow().nn_jet_index;
                        }
                        if dist < jet_b.borrow().nn_dist {
                            jet_b.borrow_mut().nn_dist = dist;
                            jet_b.borrow_mut().nn_jet_index = jet_a.borrow().nn_jet_index;
                        }
                        jet_b_next = jet_a.borrow().next_jet.clone();
                    }
                    jet_a_next = jet_a.borrow().next_jet.clone();
                }
            }
        });

        // now we want to create the DiJ distance metric where NN is J
        // this is easy since all the info is already computed prev
        // just loop over all jets
        // TODO fix this goofy ProxyJet shenanigans
        let mut di_j: Vec<f64> = vec![0.0; n];
        di_j.iter_mut().enumerate().for_each(|(i, jet)| {
            *jet = <Rc<RefCell<TiledJet>> as ProxyJet>::_bj_dij(&bj_jets[i], &bj_jets);
        });

        //start recombination loop
        for i in 0..n {
            let end_idx = n - 1 - i;
            let mut min_idx = 0;
            let mut min_dij = di_j[0];
            for (k, jet) in di_j.iter().enumerate().take(end_idx + 1).skip(1) {
                if *jet < min_dij {
                    min_dij = *jet;
                    min_idx = k;
                }
            }

            let jet_a = &bj_jets[min_idx]; //might need to clone tbd
            let mut jet_b_idx = jet_a.nn_jet_index();

            // normalize
            let min_dij = min_dij * self._invr2;

            // keep old pointer info after removing from tiles
            let mut old_b_jet: Option<Rc<RefCell<TiledJet>>> = None;

            //check if jet_a has neighbor
            match jet_b_idx {
                Some(mut _jet_b_idx) => {
                    if min_idx < _jet_b_idx {
                        std::mem::swap(&mut min_idx, &mut _jet_b_idx);
                    }
                    // nn is the new index we want to store the recombined index at
                    // modify jetB for BJ with new info

                    let new_bj = self.do_jet_jet_recombination_step::<TiledJet>(
                        bj_jets[min_idx].particle_index(),
                        bj_jets[_jet_b_idx].particle_index(),
                        min_dij,
                    );

                    // do jj now remove jetA from tiles
                    // and remove old jetB
                    // and add newly merged jetB
                    self.bj_remove_from_tiles(&bj_jets[min_idx]);
                    // no longer need this line because we call tj_set_jetinfo
                    // bj_jets[_jet_b_idx] = new_bj;
                    // apperantly need oldJetB pointer needs to be stored?
                    old_b_jet = Some(bj_jets[_jet_b_idx].clone());
                    self.bj_remove_from_tiles(&bj_jets[_jet_b_idx]);

                    self._tj_set_jetinfo(Some(new_bj), bj_jets[_jet_b_idx].particle_index())
                        .particle_index();

                    // TODO: check if removing this line doesnt do anything
                    jet_b_idx = Some(_jet_b_idx);
                }
                None => {
                    self.do_jet_beam_recombination_step(bj_jets[min_idx].particle_index(), min_dij);

                    //remove jet from tile
                    self.bj_remove_from_tiles(&bj_jets[min_idx])
                }
            };

            //some tile stuff check in between to find out which neighbors we need to update

            let mut n_near_tiles: usize = 0;
            // add neighbors to one search vector "tile_union"
            self._add_neighbors_to_tile_union(
                bj_jets[min_idx].borrow().tile_index,
                &mut tile_union,
                &mut n_near_tiles,
            );

            if let Some(jet_b_idx) = jet_b_idx {
                // if we add other neighbors then add flag
                let mut need_sort = false;
                // if same tile then ignore
                if bj_jets[jet_b_idx].borrow().tile_index != bj_jets[min_idx].borrow().tile_index {
                    need_sort = true;
                    self._add_neighbors_to_tile_union(
                        bj_jets[jet_b_idx].borrow().tile_index,
                        &mut tile_union,
                        &mut n_near_tiles,
                    );
                }

                //check if old jet needs to be added too
                if let Some(old_b_jet) = old_b_jet
                    && old_b_jet.borrow().tile_index != bj_jets[jet_b_idx].borrow().tile_index
                    && old_b_jet.borrow().tile_index != bj_jets[min_idx].borrow().tile_index
                {
                    need_sort = true;
                    self._add_neighbors_to_tile_union(
                        old_b_jet.borrow().tile_index,
                        &mut tile_union,
                        &mut n_near_tiles,
                    );
                }

                if need_sort {
                    // sort tiles up to n_near tiles just take a slice
                    // inplace?
                    let _ = &mut tile_union[0..n_near_tiles].sort();

                    //take sorted list and now compress repeated tiles
                    // by doing On^2 search for duplicates
                    let mut compressed_idx = 1;
                    for idx in 1..n_near_tiles {
                        if tile_union[idx] != tile_union[compressed_idx - 1] {
                            tile_union[compressed_idx] = tile_union[i];
                            compressed_idx += 1;
                        }
                    }
                    n_near_tiles = compressed_idx;
                }

                //here we need old jet pointer info
            }

            // now we update the nn neighbor info and move old jet

            //moves recombined min_idx with potential new NN jet at end_idx
            //copy tail values to jetA

            if min_idx != end_idx {
                bj_jets.swap(min_idx, end_idx);
                di_j.swap(min_idx, end_idx);

                //check entire list and update NN of references to new jet
                //very slow O(n) each time, compared to implicit pointer copying in C++
                for jet in bj_jets.iter_mut().take(end_idx) {
                    if jet.nn_jet_index() == Some(end_idx) {
                        jet.set_nn_jet(Some(min_idx));
                    }
                }
            }

            // now update new jetB NN by using tile union
            for &tile_union_neighbor in tile_union.iter().take(n_near_tiles) {
                let tile = &self.tiles_struct._tiles[tile_union[tile_union_neighbor]];

                let mut jet_i_next = tile.head.clone();
                while let Some(jet_i) = jet_i_next {
                    // if jet_i had NN of jetA or jetB need to recalculate
                    //let mut jet_b_next = tile.head.clone();
                    if jet_i.borrow().nn_jet_index == bj_jets[min_idx].borrow().nn_jet_index
                        && let Some(jet_b_idx) = jet_b_idx
                        && jet_i.borrow().nn_jet_index == bj_jets[jet_b_idx].borrow().nn_jet_index
                    {
                        jet_i.borrow_mut().nn_dist = self.r2;
                        jet_i.borrow_mut().nn_jet_index = None;
                        // now find new NN in all tiles
                        for near_tile_idx in tile.begin_tiles {
                            // now go through bj over all jets in each tile
                            let mut jet_j_next =
                                self.tiles_struct._tiles[near_tile_idx].head.clone();
                            while let Some(jet_j) = jet_j_next {
                                let dist = TiledJet::_bj_dist(&jet_i.borrow(), &jet_j.borrow());
                                // let dist = <Rc<RefCell<TiledJet>> as ProxyJet>::_bj_dist(
                                //     &jet_i, &jet_j,
                                // );
                                if dist < jet_i.borrow().nn_dist && !Rc::ptr_eq(&jet_i, &jet_j) {
                                    jet_i.borrow_mut().nn_dist = dist;
                                    // TODO: face same problem with jet_j since we dont have index
                                    // for tiled jet it may be wiser to just set nn_jet with RC now
                                    // instead of jet_index
                                    jet_i.borrow_mut().nn_jet_index =
                                        Some(jet_j.borrow().bj_jet_index);
                                }
                                jet_j_next = jet_j.borrow().next_jet.clone();
                            }
                        }
                    }
                    // if JJ recomb check if new jet is closer than jeti's current NN
                    if let Some(jet_b_idx) = jet_b_idx {
                        let dist =
                            TiledJet::_bj_dist(&jet_i.borrow(), &bj_jets[jet_b_idx].borrow());
                        // let dist = <Rc<RefCell<TiledJet>> as ProxyJet>::_bj_dist(
                        //     &jet_i,
                        //     &bj_jets[jet_b_idx],
                        // );
                        if dist < jet_i.borrow().nn_dist && !Rc::ptr_eq(&jet_i, &bj_jets[jet_b_idx])
                        {
                            jet_i.borrow_mut().nn_dist = dist;
                            jet_i.borrow_mut().nn_jet_index = Some(jet_b_idx);
                            // how to get index of jet_i for di_j if jet_i in this loop has no index????
                            // TODO figure out jet_i index for
                            di_j[jet_i.borrow().bj_jet_index] =
                                <Rc<RefCell<TiledJet>> as ProxyJet>::_bj_dij(&jet_i, &bj_jets); // update diJ...
                        }
                        if dist < bj_jets[jet_b_idx].borrow().nn_dist
                            && !Rc::ptr_eq(&jet_i, &bj_jets[jet_b_idx])
                        {
                            bj_jets[jet_b_idx].borrow_mut().nn_dist = dist;
                            bj_jets[jet_b_idx].borrow_mut().nn_jet_index = Some(jet_b_idx);
                            //do not need dij calc here
                        }
                    }

                    // cloning here is proper since its just an Rc but also slow
                    jet_i_next = jet_i.borrow().next_jet.clone();
                }
            }

            // update diJ for our new jetB
            if let Some(jet_b_idx) = jet_b_idx {
                di_j[jet_b_idx] =
                    <Rc<RefCell<TiledJet>> as ProxyJet>::_bj_dij(&bj_jets[jet_b_idx], &bj_jets);
            }

            // update all pointers of old tail tiles
            for near_tiled_idx in
                self.tiles_struct._tiles[bj_jets[end_idx].borrow().tile_index].begin_tiles
            {
                let near_tile = &self.tiles_struct._tiles[near_tiled_idx];
                let mut jet_j_next = near_tile.head.clone();
                while let Some(jet_j) = jet_j_next {
                    if jet_j.borrow().nn_jet_index == Some(end_idx) {
                        jet_j.borrow_mut().nn_jet_index = Some(min_idx);
                    }
                    jet_j_next = jet_j.borrow().next_jet.clone();
                }
            }
        }
    }

    #[inline]
    #[cfg(feature = "simd")]
    fn set_nn_bj_vec_simd_and_dij<J: ProxyJet>(&self, bj_jets: &mut Vec<J>, di_j: &mut Vec<f64>) {
        // can obtain simd vectors from slices, just need to handle remainder
        // for non SIMD vectors we need to init one whole array to vectorize
        // for remainder we can just ignore and run regular scalar algorithm

        // however if jet does not have NN we must set to f64::max for getting min
        // thus we need nn_array, nn_dist array
        // range to consider
        // rap_array, phi_array, R2 as splat (constant)
        // handle the dist function with these array

        let n = bj_jets.len();
        let mut rap_array = Vec::with_capacity(n);
        let mut phi_array = Vec::with_capacity(n);
        let mut nn_array = Vec::with_capacity(n);
        let mut nn_dist_array = Vec::with_capacity(n);
        let mut kt2_array = Vec::with_capacity(n);

        // push all vectors now
        bj_jets.iter().enumerate().for_each(|(i, jet)| {
            rap_array.push(jet.eta());
            phi_array.push(jet.phi());
            nn_array.push(jet.nn_jet_index().unwrap_or_else(|| i) as u64);
            nn_dist_array.push(jet.nn_dist());
            kt2_array.push(jet.kt2());
        });

        // now we can acc call SIMD func
        //
        // define a const for lanes

        const LANES: usize = 16;
        let n = bj_jets.len();
        for i in 1..n {
            ClusterSequence::set_nn_crosscheck_simd::<LANES>(
                i,
                0,     // from = 0
                i - 1, // to s= i-1 inclusive
                &rap_array,
                &phi_array,
                self.r2,
                &mut nn_dist_array,
                &mut nn_array,
            );
            // // we still need to handle remainder
            // // j does not return i
            // for remainder_idx in j..=i {
            //     self.bj_set_nn_crosscheck(&mut bj_jets[0..=remainder_idx]);
            // }
        }

        ClusterSequence::_dij_simd::<LANES>(&kt2_array, &nn_array, &nn_dist_array, di_j);

        //now go through all bj jets and set nn and nndist
        bj_jets.iter_mut().enumerate().for_each(|(i, jet)| {
            jet.set_nn_dist(nn_dist_array[i]);
            jet.set_nn_jet(if nn_array[i] as usize != i {
                Some(nn_array[i] as usize)
            } else {
                None
            });
        });
    }

    #[inline]
    #[cfg(feature = "simd")]
    pub fn set_nn_crosscheck_simd<const LANES: usize>(
        i: usize,
        start: usize,
        end: usize, // inclusive; if from > to, do nothing
        rap_arr: &[f64],
        phi_arr: &[f64],
        r2: f64,
        nndist_arr: &mut [f64],
        nn_arr: &mut [u64],
    ) {
        let rap_i = rap_arr[i];

        let phi_i = phi_arr[i];

        let mut j = start;

        // SIMD state: best candidate for i across all processed j's
        let mut best_d2 = Simd::<f64, LANES>::splat(r2);
        let mut best_j = Simd::<u64, LANES>::splat(i as u64);

        // constant lane offsets [0,1,2,...]
        let lane_offsets: Simd<u64, LANES> = Simd::from_array(std::array::from_fn(|k| k as u64));

        while j + LANES <= (end + 1) {
            let rap_j = Simd::<f64, LANES>::from_slice(&rap_arr[j..j + LANES]);
            let phi_j = Simd::<f64, LANES>::from_slice(&phi_arr[j..j + LANES]);
            let dist = ClusterSequence::bj_dist_simd(rap_i, phi_i, rap_j, phi_j);
            // 1) Update best candidate for i (vector select)
            let idx = Simd::<u64, LANES>::splat(j as u64) + lane_offsets;
            let better_for_i = dist.simd_lt(best_d2);
            best_d2 = better_for_i.select(dist, best_d2);
            best_j = better_for_i.select(idx, best_j);

            // 2) Cross-check update for [j..j+LANES): nndist/nn, via masked blend
            //    Because this region is contiguous, we can do SIMD loads/stores.
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

        //handle the tail/remainder of non simd aligned lanes regularly
        while j <= end {
            let drap = rap_i - rap_arr[j];
            let mut dphi = (phi_i - phi_arr[j]).abs();
            if dphi > PI {
                dphi = 2.0 * PI - dphi
            }
            let d2 = drap.mul_add(drap, dphi * dphi);

            // update i's best (scalar for tail)
            if d2 < nndist_min {
                nndist_min = d2;
                nn_min = j as u64;
            }

            // cross-check update (scalar for tail)
            if d2 < nndist_arr[j] {
                nndist_arr[j] = d2;
                nn_arr[j] = i as u64;
            }

            j += 1;
        }

        //Final horizontal reduction once
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

    #[inline]
    #[cfg(feature = "simd")]
    fn bj_dist_simd<const LANES: usize>(
        rap_i: f64,
        phi_i: f64,
        rap_j: Simd<f64, LANES>,
        phi_j: Simd<f64, LANES>,
    ) -> Simd<f64, LANES> {
        //splat to match lane size
        let rap_i = Simd::splat(rap_i);
        let phi_i = Simd::splat(phi_i);

        let drap = rap_i - rap_j;

        let mut dphi = (phi_i - phi_j).abs();
        let pi = Simd::splat(PI);
        let two_pi = Simd::splat(2.0 * PI);

        let mask = dphi.simd_gt(pi);
        dphi = mask.select(two_pi - dphi, dphi);

        dphi.mul_add(dphi, drap * drap)
    }

    #[inline]
    #[cfg(feature = "simd")]
    fn _dij_simd<const LANES: usize>(
        kt2_arr: &[f64],  // All kt2 values
        nn_idx: &[u64],   // Neighbor indices
        nn_dist: &[f64],  // Neighbor distances
        di_j: &mut [f64], // Output
    ) {
        let n = kt2_arr.len();
        let mut j = 0;

        while j + LANES <= n {
            // 1. Load data from slices into registers
            let kt2_i = Simd::<f64, LANES>::from_slice(&kt2_arr[j..j + LANES]);
            let dist_i = Simd::<f64, LANES>::from_slice(&nn_dist[j..j + LANES]);
            let indices = Simd::<u64, LANES>::from_slice(&nn_idx[j..j + LANES]);

            let indices_usize = indices.cast::<usize>();

            // 2. The Gather: kt2_b = jets[index].kt2()
            // We use gather_or to look up kt2_arr[indices[0]], kt2_arr[indices[1]], etc.
            // If an index is out of bounds (shouldn't happen), it defaults to kt2_i.
            let kt2_b = Simd::<f64, LANES>::gather_or(kt2_arr, indices_usize, kt2_i);

            // 3. The Math: min(kt2_i, kt2_b) * dist_i
            // Replacing 'if kt2_b < kt2_i' with a SIMD min
            let kt2_min = kt2_i.simd_min(kt2_b);
            let res = kt2_min * dist_i;

            // 4. Store result back to slice
            res.copy_to_slice(&mut di_j[j..j + LANES]);

            j += LANES;
        }

        // Scalar tail handling
        for i in j..n {
            let idx = nn_idx[i] as usize;
            di_j[i] = nn_dist[i] * kt2_arr[i].min(kt2_arr[idx]);
        }
    }

    #[inline(always)]
    #[cfg(not(feature = "simd"))]
    fn set_nn_and_dij<J: ProxyJet>(&mut self, bj_jets: &mut Vec<J>, di_j: &mut Vec<f64>) {
        let n = bj_jets.len();
        for i in 1..n {
            self.bj_set_nn_crosscheck(&mut bj_jets[0..=i]);
        }

        di_j.iter_mut().enumerate().for_each(|(i, jet)| {
            *jet = J::_bj_dij(&bj_jets[i], &bj_jets);
        });
    }
    #[cfg(feature = "simd")]
    fn set_nn_and_dij<J: ProxyJet>(&mut self, bj_jets: &mut Vec<J>, di_j: &mut Vec<f64>) {
        self.set_nn_crosscheck_simd(bj_jets, di_j);
    }

    // TODO: change from briefjets to generic ProxyJet class
    pub fn simple_n2_cluster<J: ProxyJet>(&mut self) {
        let n = self.particles.len();

        let mut bj_jets: Vec<J> = Vec::with_capacity(n);

        // set brief jet info
        // TODO; prolly will change when multithreading because of iterator and vec push
        self.particles.iter().enumerate().for_each(|(i, jet)| {
            let bj: J = ProxyJet::_bj_set_jetinfo(
                i,
                jet,
                self.r2,
                ClusterSequence::jet_scale_for_algorithm(
                    jet,
                    &self.jetdef.algorithm,
                    self.jetdef.extra_param,
                ),
            );
            bj_jets.push(bj);
        });

        // only add cfg_if when we acc implement simd feature
        // cfg_if! {
        //         if #[cfg(crate = "f64_SIMD")] {
        //             set_nn_info_simd(&mut bj_jets[0..=i], i)
        //         } else {
        //             for i in 1..n {
        //                 self.bj_set_nn_crosscheck(&mut bj_jets[0..=i]);
        //             }
        //         }
        //     }
        // // cand for simd

        // self.set_nn_bj_vec_simd(&mut bj_jets);
        //
        //self.set_nn_bj_vec_simd_and_dij(&mut bj_jets);

        // for i in 1..n {
        //     self.bj_set_nn_crosscheck(&mut bj_jets[0..=i]);
        // }

        let mut di_j: Vec<f64> = vec![0.0; n];
        self.set_nn_and_dij(&mut bj_jets, &mut di_j);

        // di_j.iter_mut().enumerate().for_each(|(i, jet)| {
        //     *jet = J::_bj_dij(&bj_jets[i], &bj_jets);
        // });

        for i in 0..n {
            //find min distance
            let end_idx = n - 1 - i;
            let mut min_idx = 0;
            let mut min_dij = di_j[0];
            for (k, jet) in di_j.iter().enumerate().take(end_idx + 1).skip(1) {
                if *jet < min_dij {
                    min_dij = *jet;
                    min_idx = k;
                }
            }

            let jet_a = &bj_jets[min_idx];
            let mut jet_b_idx = jet_a.nn_jet_index(); //THIS IS AN OPTION

            // normalize
            let min_dij = min_dij * self._invr2;

            //check if jet_a has neighbor
            match jet_b_idx {
                Some(mut _jet_b_idx) => {
                    if min_idx < _jet_b_idx {
                        std::mem::swap(&mut min_idx, &mut _jet_b_idx);
                    }
                    // nn is the new index we want to store the recombined index at
                    // modify jetB for BJ with new info

                    let new_bj = self.do_jet_jet_recombination_step(
                        bj_jets[min_idx].particle_index(),
                        bj_jets[_jet_b_idx].particle_index(),
                        min_dij,
                    );
                    bj_jets[_jet_b_idx] = new_bj;
                    jet_b_idx = Some(_jet_b_idx);
                }
                None => {
                    self.do_jet_beam_recombination_step(bj_jets[min_idx].particle_index(), min_dij);
                }
            };

            //moves recombined min_idx with potential new NN jet at end_idx
            //copy tail values to jetA
            bj_jets.swap(min_idx, end_idx);
            di_j.swap(min_idx, end_idx);

            //check entire list and update NN of references to new jet
            //very slow O(n) each time, compared to implicit pointer copying in C++
            for jet in bj_jets.iter_mut().take(end_idx) {
                if jet.nn_jet_index() == Some(end_idx) {
                    jet.set_nn_jet(Some(min_idx));
                }
            }

            // find all NN of recombined jet_b
            let tail = end_idx;
            match jet_b_idx {
                Some(_jet_b_idx) => {
                    for i in 0..tail {
                        // see if jetI had jetA or jetB as a NN if so recalc NN
                        if bj_jets[i].nn_jet_index() == Some(min_idx)
                            || bj_jets[i].nn_jet_index() == jet_b_idx
                        {
                            self.bj_set_nn_nocross(i, 0, tail, &mut bj_jets[0..tail]);
                            di_j[i] = J::_bj_dij(&bj_jets[i], &bj_jets);
                        }

                        if _jet_b_idx == i {
                            let jet_i = &mut bj_jets[i];
                            if jet_i.nn_jet_index() == Some(tail) {
                                jet_i.set_nn_jet(Some(min_idx));
                            }
                            continue;
                        }

                        let jet_ib_dist = J::_bj_dist(&bj_jets[i], &bj_jets[_jet_b_idx]);

                        if jet_ib_dist < bj_jets[i].nn_dist() {
                            bj_jets[i].set_nn_dist(jet_ib_dist);
                            bj_jets[i].set_nn_jet(jet_b_idx);
                            di_j[i] = J::_bj_dij(&bj_jets[i], &bj_jets);
                        }

                        // check if jetI is potentially jetB's NN
                        if jet_ib_dist < bj_jets[_jet_b_idx].nn_dist() {
                            bj_jets[_jet_b_idx].set_nn_dist(jet_ib_dist);
                            bj_jets[_jet_b_idx].set_nn_jet(Some(i));
                        }

                        // if jetI's NN is the new tail then relabel so it becomes jetA?
                        if bj_jets[i].nn_jet_index() == Some(tail) {
                            bj_jets[i].set_nn_jet(Some(min_idx));
                        }
                    }
                    // update new bj_dij for jetB in diJ arr
                    di_j[_jet_b_idx] = J::_bj_dij(&bj_jets[_jet_b_idx], &bj_jets);
                }
                None => {
                    //check if old JetA was NN to jet being mapped
                    for i in 0..tail {
                        if let Some(jet_i_nn_idx) = bj_jets[i].nn_jet_index() {
                            if jet_i_nn_idx == min_idx {
                                self.bj_set_nn_nocross(i, 0, tail, &mut bj_jets[0..tail]);
                                di_j[i] = J::_bj_dij(&bj_jets[i], &bj_jets);
                            }
                            //if jet has NN of tail then used to be jetA
                            if jet_i_nn_idx == end_idx {
                                bj_jets[jet_i_nn_idx].set_nn_jet(Some(min_idx));
                            }
                        }
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
                "Child {} history index for {} already set",
                self._history[parent1].child, parent1
            );
        }
        self._history[parent1].child = child_hist_index;

        if !(parent2 == BEAMJET || parent2 == INVALID || parent2 == INEXISTENT_PARENT) {
            if self._history[parent2].child != INVALID {
                panic!(
                    "Child {} history index for {} already set",
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
        jets: &mut [J],
    ) {
        let mut nn_dist = self.r2;
        let mut nn: Option<usize> = None;

        // basically splitting for loop into two parts
        // can prolly use this in for loop with a || stop condition
        if head_idx < curr_idx {
            for jet_b_idx in head_idx..curr_idx {
                let dist = J::_bj_dist(&jets[curr_idx], &jets[jet_b_idx]);
                if dist < nn_dist {
                    nn_dist = dist;
                    nn = Some(jet_b_idx);
                }
            }
        }

        if tail_idx > curr_idx {
            for jet_b_idx in curr_idx + 1..tail_idx {
                let dist = J::_bj_dist(&jets[curr_idx], &jets[jet_b_idx]);
                if dist < nn_dist {
                    nn_dist = dist;
                    nn = Some(jet_b_idx);
                }
            }
        }

        jets[curr_idx].set_nn_jet(nn);
        jets[curr_idx].set_nn_dist(nn_dist);
    }

    #[inline]
    fn bj_set_nn_crosscheck<J: ProxyJet>(&self, jets: &mut [J]) {
        let mut nn_dist = self.r2;
        let mut nn: Option<usize> = None;
        let n = jets.len() - 1;
        for i in 0..n {
            let dist = J::_bj_dist(&jets[i], &jets[n]);
            if dist < nn_dist {
                nn_dist = dist;
                nn = Some(i);
            }
            if dist < jets[i].nn_dist() {
                jets[i].set_nn_dist(dist);
                jets[i].set_nn_jet(Some(n));
            }
        }
        jets[n].set_nn_dist(nn_dist);
        jets[n].set_nn_jet(nn);
    }
}
