use log::warn;
use std::panic;

use crate::constants::PI;
use crate::proxy_jet::BriefJet;
use crate::proxy_jet::ProxyJet;
use crate::proxy_jet::Tile;
//use crate::proxy_jet::TiledJet;
use crate::pseudo_jet::PseudoJet;

// TODO: add automatic strategies for clustering like BEST: N2Plain, N2Tiling
pub struct ClusterSequence {
    pub particles: Vec<PseudoJet>,
    pub jetdef: JetDefinition,
    _invr2: f64,
    r2: f64,
    _history: Vec<HistoryElement>,
    init_n: usize,
    //TODO: find if we can move tiled jets related vars out of here
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
}

// TODO: implement other algorithms
#[derive(PartialEq, Eq)]
pub enum Algorithm {
    AntiKt,
    Kt,
    Cambridge,
    Genkt,
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
            //TODO: remvoe all this goofy ahh inits
            _tiles: Vec::new(),
            _tiles_eta_min: 0.0,
            _tiles_eta_max: 0.0,
            _tile_size_eta: 0.0,
            _tile_size_phi: 0.0,
            _n_tiles_phi: 0,
            _tiles_ieta_min: 0,
            _tiles_ieta_max: 0,
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
            Algorithm::Genkt => {
                let p = extra_param.unwrap_or(0.0);
                // TODO: figure out if there's a better way to handle this
                if p <= 0.0 && kt2 < 1e-300 {
                    kt2 = 1e-300;
                }
                kt2.powf(p)
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

        //TODO: when implementing other strategies have to handle here
        // since we only implemented kt and N2 this func is pretty naive

        match self.jetdef.strategy {
            Strategy::Best => {
                //TODO: implement best strategy
            }
            Strategy::N2Plain => {
                self.simple_n2_cluster();
            }
            Strategy::N2Tiling => {
                self.tiled_n2_cluster();
                //TODO: implement N2Tiling strategy
            } // _ => {
              //     panic!("Unsupported strategy");
              // }
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
            Algorithm::Genkt | Algorithm::AntiKt => self
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
            && *algorithm != Algorithm::Genkt
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

            //TODO: check if all this casting is best
            // rap will be a float so cast to i64 for binning
            let ibin = rap as i64 + nrap as i64;
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
        for i in 0..nbins {
            cum_lo += counts[i] as f64;
            if cum_lo >= allowed_max_cumul {
                let y = (i - nrap) as f64;
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
                let y = (i - nrap + 1) as f64;
                if y < max_rap {
                    max_rap = y;
                    ubin = i;
                }
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
            for i in ibin_lo + 1..ibin_hi {
                cumul_2 += counts[i] as f64 * counts[i] as f64;
            }
        }

        (min_rap, max_rap, cumul_2)
    }

    //----------------------------------------------------------------------
    /// Set up the tiles:
    ///  - decide the range in eta
    ///  - allocate the tiles
    ///  - set up the cross-referencing info between tiles
    ///
    /// The neighbourhood of a tile is set up as follows
    ///
    /// 	      LRR
    ///           LXR
    ///           LLR
    ///
    /// such that tiles is an array containing XLLLLRRRR with pointers
    ///                                         |   \ RH_tiles
    ///                                         \ surrounding_tiles
    ///
    /// with appropriate precautions when close to the edge of the tiled
    /// region.
    ///
    //setup correct tile size and point tiles to their neighbors
    pub fn initialize_tiles(&mut self) {
        // want to bound this for very low DeltaR to avoid memory blow ups
        let default_tile_size: f64 = self.jetdef.r.max(0.1);
        self._tile_size_eta = default_tile_size;

        // for phi we need to check spacing against 2pi
        // when phi is <3 no tiling is done so min = 3
        self._n_tiles_phi = ((PI * 2.0 / default_tile_size).floor() as isize).max(3);
        self._tile_size_phi = PI * 2.0 / self._n_tiles_phi as f64;

        // find the min and max rap/eta that we should be using for this analysis
        let (min_rap, max_rap, _cumul_2) = self.min_max_rap();

        //find min/max values for these tiles and figure out why _tiles_eta_min isnt just min_rap?
        self._tiles_ieta_min = (min_rap / self._tile_size_eta).floor() as isize;
        self._tiles_ieta_max = (max_rap / self._tile_size_eta).floor() as isize;
        self._tiles_eta_min = self._tiles_ieta_min as f64 * self._tile_size_eta;
        self._tiles_eta_max = self._tiles_ieta_max as f64 * self._tile_size_eta;

        // allocate the vector for the size of tiles before pushing
        // TODO: figure out how to access internal _tiles struct with object
        // (should it just be defined as TiledJet object and no trait for PseudoJet)
        // TODO: C++ code does not need to init vectors to get pointers could be perf loss
        let mut tiles: Vec<Tile> = Vec::new();
        tiles.resize(
            (self._tiles_ieta_max - self._tiles_ieta_min + 1) as usize * self._n_tiles_phi as usize,
            Tile::new(),
        );

        // now link all these tiles together
        for ieta in self._tiles_ieta_min..self._tiles_ieta_max {
            for iphi in 0..self._n_tiles_phi {
                let mut tile_idx = 0;
                let tile: &mut Tile = &mut tiles[self._tile_index(ieta, iphi)];
                // first tile has no HEAD
                tile.head = Option::None;

                //begin_tiles ppoints to neighboring tiles including itself?
                tile.begin_tiles[tile_idx] = self._tile_index(ieta, iphi);

                // the surrounding tiles excludes self
                // can just use slices in rust
                tile_idx += 1;
                //tile.surrounding_tiles = &tile.begin_tiles[tile_idx..];

                // after first loop of ieta we can now check prev column
                // now left right middle
                if ieta > self._tiles_ieta_min {
                    for idphi in [-1isize, 0, 1] {
                        tile_idx += 1;
                        tile.begin_tiles[tile_idx] = self._tile_index(ieta - 1, iphi + idphi);
                    }
                }

                // note: phi doesnt need bound checks since it mod 2pi
                // LHS includes under curr elem
                tile_idx += 1;
                tile.begin_tiles[tile_idx] = self._tile_index(ieta, iphi - 1);

                // first RHS is above curr elem
                tile_idx += 1;
                tile.begin_tiles[tile_idx] = self._tile_index(ieta, iphi + 1);

                // only set last R if we are not at max
                if ieta < self._tiles_ieta_max {
                    for idphi in [-1isize, 0, 1] {
                        tile_idx += 1;
                        tile.begin_tiles[tile_idx] = self._tile_index(ieta + 1, iphi + idphi);
                    }
                }

                // TODO: since i think this is unitialized in prev config , aka at correct points
                // check if this does not have runtime errors?
                tile.surrounding_tiles = 1..5;
                tile.rh_tiles = 5..9;
            }
        }
    }

    // get index even if iphi is negative since wraparound
    fn _tile_index(&self, ieta: isize, iphi: isize) -> usize {
        // use mod to get wraparound behavior but -1 mod n is == -1 so add by n
        ((ieta - self._tiles_ieta_min) * self._n_tiles_phi
            + ((iphi + self._n_tiles_phi) % self._n_tiles_phi)) as usize
    }

    pub fn tiled_n2_cluster(&mut self) {}

    pub fn simple_n2_cluster(&mut self) {
        let n = self.particles.len();

        let mut bj_jets: Vec<BriefJet> = Vec::with_capacity(n);

        // set brief jet info
        // TODO; prolly will change when multithreading because of iterator and vec push
        self.particles.iter().enumerate().for_each(|(i, jet)| {
            let kt2 = ClusterSequence::jet_scale_for_algorithm(
                jet,
                &self.jetdef.algorithm,
                self.jetdef.extra_param,
            );
            let bj = BriefJet {
                eta: *(jet.rap()),
                phi: *(jet.phi()),
                kt2,
                _jets_index: i,
                nn_dist: self.r2,
                nn_jet_index: None,
            };
            bj_jets.push(bj);
        });

        for i in 1..n {
            self.bj_set_nn_crosscheck(&mut bj_jets[0..=i]);
        }

        let mut di_j: Vec<f64> = vec![0.0; n];

        di_j.iter_mut().enumerate().for_each(|(i, jet)| {
            *jet = BriefJet::_bj_dij(&bj_jets[i], &bj_jets);
        });

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
                        bj_jets[min_idx]._jets_index,
                        bj_jets[_jet_b_idx]._jets_index,
                        min_dij,
                    );
                    bj_jets[_jet_b_idx] = new_bj;
                    jet_b_idx = Some(_jet_b_idx);
                }
                None => {
                    self.do_jet_beam_recombination_step(bj_jets[min_idx]._jets_index, min_dij);
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
                            di_j[i] = BriefJet::_bj_dij(&bj_jets[i], &bj_jets);
                        }

                        if _jet_b_idx == i {
                            let jet_i = &mut bj_jets[i];
                            if jet_i.nn_jet_index() == Some(tail) {
                                jet_i.set_nn_jet(Some(min_idx));
                            }
                            continue;
                        }

                        let jet_ib_dist = BriefJet::_bj_dist(&bj_jets[i], &bj_jets[_jet_b_idx]);

                        if jet_ib_dist < bj_jets[i].nn_dist() {
                            bj_jets[i].set_nn_dist(jet_ib_dist);
                            bj_jets[i].set_nn_jet(jet_b_idx);
                            di_j[i] = BriefJet::_bj_dij(&bj_jets[i], &bj_jets);
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
                    di_j[_jet_b_idx] = BriefJet::_bj_dij(&bj_jets[_jet_b_idx], &bj_jets);
                }
                None => {
                    //check if old JetA was NN to jet being mapped
                    for i in 0..tail {
                        if let Some(jet_i_nn_idx) = bj_jets[i].nn_jet_index() {
                            if jet_i_nn_idx == min_idx {
                                self.bj_set_nn_nocross(i, 0, tail, &mut bj_jets[0..tail]);
                                di_j[i] = BriefJet::_bj_dij(&bj_jets[i], &bj_jets);
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

    fn do_jet_jet_recombination_step(
        &mut self,
        jet_a_idx: usize,
        jet_b_idx: usize,
        min_dij: f64,
    ) -> BriefJet {
        let new_jet = self
            .jetdef
            .recombine(&self.particles[jet_a_idx], &self.particles[jet_b_idx]);
        let kt2 = ClusterSequence::jet_scale_for_algorithm(
            &new_jet,
            &self.jetdef.algorithm,
            self.jetdef.extra_param,
        );
        let new_bj = BriefJet {
            eta: *(new_jet.rap()),
            phi: *(new_jet.phi()),
            kt2,
            _jets_index: self.particles.len(),
            nn_dist: self.r2,
            nn_jet_index: None,
        };

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
                let dist = BriefJet::_bj_dist(&jets[curr_idx], &jets[jet_b_idx]);
                if dist < nn_dist {
                    nn_dist = dist;
                    nn = Some(jet_b_idx);
                }
            }
        }

        if tail_idx > curr_idx {
            for jet_b_idx in curr_idx + 1..tail_idx {
                let dist = BriefJet::_bj_dist(&jets[curr_idx], &jets[jet_b_idx]);
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
            let dist = BriefJet::_bj_dist(&jets[i], &jets[n]);
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
