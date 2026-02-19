use std::panic;
use std::usize;

use crate::constants::PI;
use crate::proxy_jet::BriefJet;
use crate::proxy_jet::JetType;
use crate::proxy_jet::ProxyJet;
use crate::pseudo_jet::PseudoJet;

// TODO: add automatic strategies for clustering like BEST: N2Plain, N2Tiling
pub struct ClusterSequence {
    pub particles: Vec<PseudoJet>,
    pub jetdef: JetDefinition,
    _invr2: f64,
    r2: f64,
    _history: Vec<HistoryElement>,
    init_n: usize,
}

pub struct JetDefinition {
    pub algorithm: Algorithm,
    pub r: f64,
    pub scheme: RecombinationScheme,
    pub strategy: Strategy,
}

pub enum Strategy {
    Best,
    N2Plain,
    N2Tiling,
}

// TODO: implement other algorithms
pub enum Algorithm {
    AntiKt,
    Kt,
    Cambridge,
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
    ) -> Self {
        JetDefinition {
            algorithm,
            r,
            scheme,
            strategy,
        }
    }

    //TODO: should add reset function to make it faster
    pub fn recombine(&self, jet_a: &PseudoJet, jet_b: &PseudoJet) -> PseudoJet {
        match self.scheme {
            EScheme => {
                let jet_ab = jet_b + jet_a; //keeps mem in jet_b
                jet_ab
            }
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

// TODO: rlly need to rename this into smth more "representative" (like jettype lololol)
// normally would use neg values but usize indexing is unsigned
#[repr(usize)]
pub enum JetErrors {
    Invalid = usize::MAX,
    BeamJet = usize::MAX - 1,
    InexistentParent = usize::MAX - 2,
}

#[allow(dead_code)]
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
        };

        clust_seq.fill_init_history();

        clust_seq
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
                //TODO: implement N2Tiling strategy
            } // _ => {
              //     panic!("Unsupported strategy");
              // }
        }
    }

    fn fill_init_history(&mut self) {
        //get total energy of the event
        let mut tot_e = 0.0;
        //map for all _jets
        self.particles.iter_mut().enumerate().for_each(|(i, jet)| {
            self._history.push(HistoryElement {
                parent1: JetErrors::InexistentParent as usize,
                parent2: JetErrors::InexistentParent as usize,
                child: JetErrors::Invalid as usize,
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

    pub fn inclusive_jets(&self, pmin: f64) -> Vec<PseudoJet> {
        let dcut = pmin * pmin;
        match self.jetdef.algorithm {
            Algorithm::Kt => {
                // loop through history in reverse
                let jets = self
                    ._history
                    .iter()
                    .rev()
                    .take_while(|hist_elem| hist_elem.curr_max_d_ij >= dcut)
                    .filter_map(|hist_elem| {
                        if hist_elem.parent2 == (JetErrors::BeamJet as usize)
                            && hist_elem.d_ij >= dcut
                        {
                            Some(self.particles[self._history[hist_elem.parent1].jet_index])
                        } else {
                            None
                        }
                    })
                    .collect();
                jets
            }
            Algorithm::AntiKt => self
                ._history
                .iter()
                .rev()
                .filter_map(|hist_elem| {
                    if hist_elem.parent2 == (JetErrors::BeamJet as usize) {
                        let parent1 = hist_elem.parent1;
                        let jet = self.particles[self._history[parent1].jet_index];
                        if jet.perp2() >= dcut { Some(jet) } else { None }
                    } else {
                        None
                    }
                })
                .collect(),
            _ => panic!("Unsupported algorithm"),
        }
    }

    fn split_bj_mut_slices<'a, J: ProxyJet>(
        &self,
        i_idx: usize,
        j_idx: usize,
        jets: &'a mut [J],
    ) -> (&'a mut J, &'a mut J) {
        let split_idx = std::cmp::min(i_idx, j_idx);

        let left_bj_jet: &mut J;
        let right_bj_jet: &mut J;

        // we are splitting by i and i< jet_b_idx
        // TODO: messy split_at mut def more generic way to do this since also
        // needed in inital nn_cross_checlk
        if split_idx == i_idx {
            let (left_bj_jets, right_bj_jets) = jets.split_at_mut(split_idx + 1);
            left_bj_jet = &mut left_bj_jets[split_idx];
            right_bj_jet = &mut right_bj_jets[j_idx - split_idx - 1];
        } else {
            let (left_bj_jets, right_bj_jets) = jets.split_at_mut(split_idx + 1);
            left_bj_jet = &mut left_bj_jets[split_idx];
            right_bj_jet = &mut right_bj_jets[i_idx - split_idx - 1];
        };

        (left_bj_jet, right_bj_jet)
    }

    // main driver loop
    pub fn simple_n2_cluster(&mut self) {
        let n = self.particles.len();

        let mut bj_jets: Vec<BriefJet> = Vec::with_capacity(n);

        // set brief jet info
        // TODO; prolly will change when multithreading because of iterator and vec push
        self.particles.iter_mut().enumerate().for_each(|(i, jet)| {
            let bj = BriefJet {
                eta: jet.rap(),
                phi: jet.phi(),
                kt2: if jet.kt2() > 1e-300 {
                    1. / jet.kt2()
                } else {
                    1e300
                }, //TODO: implement jet_scale_for_algorithm
                _jets_index: i,
                nn_dist: self.r2,
                nn_jet_index: None,
            };
            bj_jets.push(bj);
        });

        for i in 1..n {
            self.bj_set_nn_crosscheck(&mut bj_jets[0..=i]);
        }

        // calc distances
        let mut di_j: Vec<f64> = vec![0.0; n];

        di_j.iter_mut().enumerate().for_each(|(i, jet)| {
            *jet = ClusterSequence::_bj_dij(&bj_jets[i], &bj_jets);
        });

        // run recombination
        // while head neq tail
        // find minimum of the diJ
        // do the recombination between neighboring jets with minimum distance
        // update bookkeeping
        // update tail and head pointers
        for i in 0..n {
            //find min distance
            let end_idx = n - 1 - i;
            let (mut min_idx, min_dij) = di_j[0..=end_idx]
                .iter()
                .enumerate()
                .min_by(|a, b| a.1.total_cmp(&b.1))
                .unwrap();

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
            bj_jets[min_idx] = bj_jets[end_idx].clone();
            // now copy over tail values to the min distances loc also
            di_j[min_idx] = di_j[end_idx];

            //check entire list and update NN of references to new jet
            //very slow O(n) each time, compared to implicit pointer copying in C++
            for i in 0..bj_jets.len() {
                if bj_jets[i].nn_jet_index() == Some(end_idx) {
                    bj_jets[i].set_nn_jet(Some(min_idx));
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
                            di_j[i] = ClusterSequence::_bj_dij(&bj_jets[i], &bj_jets);
                        }

                        if _jet_b_idx == i {
                            let jet_i = &mut bj_jets[i];
                            if jet_i.nn_jet_index() == Some(tail) {
                                jet_i.set_nn_jet(Some(min_idx));
                            }
                            continue;
                        }

                        let mut_jets = self.split_bj_mut_slices(i, _jet_b_idx, &mut bj_jets);
                        let jet_ib_dist = ClusterSequence::bj_dist(mut_jets.0, mut_jets.1);

                        //TODO clean up next two blocks of comparing and setting nn_dist
                        {
                            if jet_ib_dist < bj_jets[i].nn_dist() {
                                {
                                    let jet_i = &mut bj_jets[i];
                                    jet_i.set_nn_dist(jet_ib_dist);
                                    jet_i.set_nn_jet(jet_b_idx);
                                }
                                let jet_i = &bj_jets[i];
                                di_j[i] = ClusterSequence::_bj_dij(jet_i, &bj_jets);
                            }
                        }

                        {
                            let jet_b = &mut bj_jets[_jet_b_idx];
                            // check if jetI is potentially jetB's NN
                            if jet_ib_dist < jet_b.nn_dist() {
                                jet_b.set_nn_dist(jet_ib_dist);
                                jet_b.set_nn_jet(Some(i));
                            }
                        }

                        let jet_i = &mut bj_jets[i];
                        // if jetI's NN is the new tail then relabel so it becomes jetA?
                        if jet_i.nn_jet_index() == Some(tail) {
                            jet_i.set_nn_jet(Some(min_idx));
                        }
                    }
                    // update new bj_dij for jetB in diJ arr
                    di_j[_jet_b_idx] = ClusterSequence::_bj_dij(&bj_jets[_jet_b_idx], &bj_jets);
                }
                None => {
                    //slice of _jets
                    //check if old JetA was NN to jet being mapped
                    for i in 0..tail {
                        match bj_jets[i].nn_jet_index() {
                            Some(jet_i_nn_idx) => {
                                if jet_i_nn_idx == min_idx {
                                    //TODO: verify this aka setting di_j and how much of bjjets to include
                                    self.bj_set_nn_nocross(i, 0, tail, &mut bj_jets[0..tail]);
                                    di_j[i] = ClusterSequence::_bj_dij(&bj_jets[i], &bj_jets);
                                }
                                //if jet has NN of tail then used to be jetA
                                if jet_i_nn_idx == end_idx {
                                    bj_jets[jet_i_nn_idx].set_nn_jet(Some(min_idx));
                                }
                            }
                            None => {}
                        };
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
        // recombine jet_a and jet_b into new_jet
        let mut new_jet = self
            .jetdef
            .recombine(&self.particles[jet_a_idx], &self.particles[jet_b_idx]);

        // push new_jet onto particles
        self.particles.push(new_jet);

        // get it's index (new size of vec - 1)
        let new_idx = self.particles.len() - 1;

        // update history of jet_a and jet_b respectively
        let hist_a_idx = self.particles[jet_a_idx].cluster_hist_index();
        let hist_b_idx = self.particles[jet_b_idx].cluster_hist_index();

        self.particles[new_idx].set_cluster_hist_index(self._history.len());

        self.add_recomb_to_history(
            std::cmp::min(hist_a_idx, hist_b_idx),
            std::cmp::max(hist_a_idx, hist_b_idx),
            new_idx,
            min_dij,
        );

        BriefJet {
            eta: new_jet.rap(),
            phi: new_jet.phi(),
            kt2: if new_jet.kt2() > 1e-300 {
                1. / new_jet.kt2()
            } else {
                1e300
            }, //TODO: implement jet_scale_for_algorithm
            _jets_index: self.particles.len() - 1,
            nn_dist: self.r2,
            nn_jet_index: None,
        }
    }

    fn do_jet_beam_recombination_step(&mut self, jet_index: usize, d_ij: f64) {
        self.add_recomb_to_history(
            self.particles[jet_index].cluster_hist_index(),
            JetErrors::BeamJet as usize,
            JetErrors::Invalid as usize,
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
        // add the recombination to the history

        //create new history_element for this recomb step
        let element: HistoryElement = HistoryElement {
            parent1: parent1,
            parent2: parent2,
            jet_index: jet_index,
            d_ij: d_ij,
            child: JetErrors::Invalid as usize,
            curr_max_d_ij: f64::max(d_ij, self._history.last().unwrap().curr_max_d_ij), //TODO: verify options and stuff
        };

        self._history.push(element);

        //TODO: some sanity checks checking in clustersequence.cc
        //TODO: add flag for writeout combinations

        let child_hist_index: usize = self._history.len() - 1;
        if self._history[parent1].child != JetErrors::Invalid as usize {
            panic!("Child history index already set");
        }
        self._history[parent1].child = child_hist_index;

        if !(parent2 == JetErrors::BeamJet as usize
            || parent2 == JetErrors::Invalid as usize
            || parent2 == JetErrors::InexistentParent as usize)
        {
            if self._history[parent2].child != JetErrors::Invalid as usize {
                panic!("Child history index already set");
            }
            self._history[parent2].child = child_hist_index;
        }

        if jet_index != JetErrors::Invalid as usize {
            if jet_index == JetErrors::BeamJet as usize
                || jet_index == JetErrors::InexistentParent as usize
            {
                panic!("Jet_index has incorrect index (Beam jet or Inexistent parent)");
            }
            self.particles[jet_index].set_cluster_hist_index(child_hist_index);
        }
    }

    fn _bj_dij<J: ProxyJet>(jet: &J, jets: &[J]) -> f64 {
        let mut kt2 = jet.kt2();
        jet.nn_jet_index().map(|index| {
            if jets[index].kt2() < kt2 {
                kt2 = jets[index].kt2();
            }
        });
        jet.nn_dist() * kt2
    }

    // BriefJet == eta, phi ,kt2, NN_dist,
    // pointer to NN, and index of jet from array

    // tiledjet is linked list and extends BJ

    //params needed R (for filtering given by user), R^2 and R^-2

    // jetA is contigious in memory, do not need to use pointers

    //instead of templates (can have diff types of Jets) we can use generics

    //recomb loop
    // find min of diJ find min distance (simple)
    // merge min distance + NN into one jet

    // if we can use maps/iterators, safely, compiler does lots of optimizations

    // if no nearest neighbor jetB == Null then we just call a bookkeeping function to add to history

    // STATIC METHODS USING PROXYJET

    // it seems like functions like .eta() .phi() are part of generic trait?
    // maybe should change PseudoJet impl into a Trait instead
    // TODO: check if these functions in default trait without const hurt runtime
    #[inline]
    fn bj_dist<J: ProxyJet>(jet_a: &mut J, jet_b: &mut J) -> f64 {
        let dphi: f64 = PI - f64::abs(PI - f64::abs(jet_a.phi() - jet_b.phi()));
        let deta: f64 = jet_a.eta() - jet_b.eta();
        dphi * dphi + deta * deta
    }

    #[inline]
    fn _bj_dist<J: ProxyJet>(jet_a: &mut J, jet_b: &mut J) -> f64 {
        let dphi: f64 = PI - f64::abs(PI - f64::abs(jet_a.phi() - jet_b.phi()));
        let deta: f64 = jet_a.eta() - jet_b.eta();
        dphi * dphi + deta * deta
    }

    #[allow(dead_code)]
    #[inline]
    fn bj_dij<'a, J: ProxyJet>(jet: &J, array_jets: &[J]) -> f64 {
        let kt2 = jet.kt2();
        jet.nn_dist()
            * (match jet.nn_jet_index() {
                Some(jet_b_index) => {
                    let jet_b_kt2 = array_jets[jet_b_index].kt2();
                    if jet_b_kt2 < kt2 { jet_b_kt2 } else { kt2 }
                }
                None => jet.nn_dist() * kt2,
            })
    }

    // find the nearest neighbor for jet
    // set nn_dist tolerance to user entered r2 value
    // set nn to null
    // for each jet_b from head to tail in jet array:
    // 1. call bjdist to calculate the distance between jet and jet_b
    // 2. check if calculated dist is less than nearest neighbor distance and update nn_dist and nn to jet_b if true
    // 3. check also if calculated dist is less than jet_b's nearest neighbor distance (ie. if jet itself is a nearest neighbor of jet_b)
    // then update nn_dist and nn to jet_b if true
    // update jet's nearest neighbor to nn and jet's nn_dist to nn_dist

    #[inline]
    fn bj_set_nn_nocross<'b, 'a, J: ProxyJet>(
        &'b mut self,
        curr_idx: usize,
        head_idx: usize,
        tail_idx: usize,
        jets: &'a mut [J],
    ) -> () {
        let mut nn_dist = self.r2;
        let mut nn: Option<usize> = None;

        // basically splitting for loop into two parts
        // can prolly use this in for loop with a || stop condition
        if head_idx < curr_idx {
            for jet_b_idx in head_idx..curr_idx {
                let mut_jets = self.split_bj_mut_slices(curr_idx, jet_b_idx, jets);
                let dist = ClusterSequence::bj_dist(mut_jets.0, mut_jets.1);
                if dist < nn_dist {
                    nn_dist = dist;
                    nn = Some(jet_b_idx);
                }
            }
        }

        if tail_idx > curr_idx {
            for jet_b_idx in curr_idx + 1..tail_idx {
                let mut_jets = self.split_bj_mut_slices(curr_idx, jet_b_idx, jets);
                let dist = ClusterSequence::bj_dist(mut_jets.0, mut_jets.1);
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
    fn bj_set_nn_crosscheck<'b, 'a, J: ProxyJet>(&'b self, jets: &'a mut [J]) -> () {
        let mut nn_dist = self.r2;
        let mut nn: Option<usize> = None;
        let n = jets.len() - 1;
        for i in 0..n {
            let mut_jets = self.split_bj_mut_slices(n, i, jets);
            let dist = ClusterSequence::bj_dist(mut_jets.0, mut_jets.1);
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
