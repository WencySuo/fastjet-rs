use std::usize;

use crate::constants::PI;
use crate::proxy_jet::BriefJet;
use crate::proxy_jet::ProxyJet;
use crate::pseudo_jet::PseudoJet;

// TODO: add automatic strategies for clustering like BEST: N2Plain, N2Tiling
pub struct ClusterSequence {
    pub particles: Vec<PseudoJet>,
    pub jetdef: JetDefinition,
    _invr2: f64,
    r2: f64,
    _history: Vec<HistoryElement>,
}

pub struct JetDefinition {
    pub algorithm: Algorithm,
    pub r: f64,
    pub scheme: RecombinationScheme,
}

// TODO: implement other algorithms
pub enum Algorithm {
    AntiKt,
    Cambridge,
}

pub enum RecombinationScheme {
    EScheme,
    PTScheme,
}

impl JetDefinition {
    pub fn new(algorithm: Algorithm, r: f64, scheme: RecombinationScheme) -> Self {
        JetDefinition {
            algorithm,
            r,
            scheme,
        }
    }

    //TODO: should add reset function to make it faster
    pub fn recombine(&self, jet_a: &PseudoJet, jet_b: &PseudoJet) -> PseudoJet {
        match self.scheme {
            RecombinationScheme::EScheme => {
                let jet_ab = jet_b + jet_a; //keeps mem in jet_b
                jet_ab
            }
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

struct HistoryElement {
    parent1: usize,
    parent2: usize,
    child: usize,
    jet_index: usize,
    d_ij: f64,
    curr_max_d_ij: f64,
}

impl ClusterSequence {
    //TODO: verify if this is right place for constant
    pub fn new(particles: Vec<PseudoJet>, jetdef: JetDefinition) -> Self {
        let r2 = jetdef.r * jetdef.r;
        let _invr2 = 1.0 / r2;
        ClusterSequence {
            particles,
            jetdef,
            _invr2,
            r2,
            _history: Vec::new(),
        }
    }

    // main driver loop
    pub fn simple_n2_cluster(&mut self) {
        let n = self.particles.len();

        let mut bj_jets: Vec<BriefJet> = Vec::with_capacity(n);

        // set brief jet info
        // TODO; prolly will change when multithreading because of iterator and vec push
        self.particles.iter().enumerate().for_each(|(i, jet)| {
            let bj = BriefJet {
                eta: jet.rap(),
                phi: jet.phi(),
                kt2: jet.kt2(), //TODO: implement jet_scale_for_algorithm
                _jets_index: i,
                nn_dist: self.jetdef.r * self.jetdef.r,
                nn_jet_index: None,
            };
            bj_jets.push(bj);
        });

        for i in 1..n {
            self.bj_set_nn_crosscheck(&mut bj_jets[0..i]);
        }

        // calc distances
        let mut di_j: Vec<f64> = vec![0.0; n];

        di_j.iter_mut().enumerate().for_each(|(i, jet)| {
            *jet += ClusterSequence::_bj_dij(&bj_jets[i], &bj_jets);
        });

        // run recombination
        // while head neq tail
        // find minimum of the diJ
        // do the recombination between neighboring jets with minimum distance
        // update bookkeeping
        // update tail and head pointers
        for i in n..0 {
            //find min distance
            //
            let end_idx = n - 1 - i;

            // impossible to not have min deltaR
            let (min_idx, min_dij) = di_j[0..=end_idx]
                .iter()
                .enumerate()
                .min_by(|a, b| a.1.total_cmp(&b.1))
                .unwrap();

            let jet_a = &bj_jets[min_idx];
            let jet_b_idx = jet_a.nn_jet_index(); //THIS IS AN OPTION
            //let jet_b: Option<&BriefJet> = jet_b_idx.map(|idx| &bj_jets[idx]);

            // normalize
            let min_dij = min_dij * self._invr2;

            //check if jet_a has neighbor
            match jet_b_idx {
                Some(jet_b_idx) => {
                    // nn is the new index we want to store the recombined index at
                    //modify jetB for BJ with new info
                    self.do_jet_jet_recombination_step(min_idx, jet_b_idx, min_dij);
                }
                None => {
                    // jet beam recombination
                    self.do_jet_beam_recombination_step(min_idx, min_dij);
                }
            };

            //moves recombined min_idx with potential new NN jet at end_idx
            //swap jetA with tail vector
            bj_jets.swap(min_idx, end_idx);
            // now swap the min distances loc also
            di_j.swap(min_idx, end_idx);

            // find all NN of recombined jet_b
            let tail = end_idx - 1;
            match jet_b_idx {
                Some(_jet_b_idx) => {
                    for i in 0..=tail {
                        // see if jetI had jetA or jetB as a NN if so recalc NN
                        
                        if bj_jets[i].nn_jet_index() == Some(min_idx)
                            || bj_jets[i].nn_jet_index() == jet_b_idx
                        {
                            self.bj_set_nn_nocross(i, 0 , tail, &mut bj_jets[0..tail]);
                        }
                        
                        if _jet_b_idx != i {
                            return
                        }
                        
                        // fight borrow checker
                        let mid = (i + _jet_b_idx) / 2; 
                        let modulo = (i + _jet_b_idx) % 2;
                        let (left_bj_jets, right_bj_jets) = bj_jets.split_at_mut(mid+modulo);
        
                        let jet_i = &mut left_bj_jets[i];
                        let jet_b = &mut right_bj_jets[_jet_b_idx - mid - modulo];
                        
                        // check if new jetB is closer than jet I's current N and update
                        let jet_ib_dist = ClusterSequence::bj_dist(jet_i, jet_b);
                        if jet_ib_dist < jet_i.nn_dist()  {
                            jet_i.set_nn_dist(jet_ib_dist);
                            jet_i.set_nn_jet(jet_b_idx);
                        }
                        // check if jetI is potentially jetB's NN
                        if jet_ib_dist < jet_b.nn_dist() {
                            jet_b.set_nn_dist(jet_ib_dist);
                            jet_b.set_nn_jet(Some(i));
                        }
                        // if jetI's NN is the new tail then relabel so it becomes jetA?
                        if jet_i.nn_jet_index() == Some(tail) {
                            jet_i.set_nn_jet(Some(min_idx));
                        }
                    };
                    // update new bj_dij for jetB in diJ arr
                    di_j[jet_b_idx.unwrap()] = ClusterSequence::_bj_dij(&bj_jets[_jet_b_idx], &bj_jets);
                }
                None => {
                    //slice of _jets
                    //check if old JetA was NN to jet being mapped
                    for i in 0..=tail {
                        match bj_jets[i].nn_jet_index() {
                            Some(jet_i_idx) => {
                                if jet_i_idx == min_idx {
                                    //TODO: verify this aka setting di_j and how much of bjjets to include
                                    self.bj_set_nn_nocross(
                                        jet_i_idx,
                                        0,
                                        tail,
                                        &mut bj_jets[0..tail],
                                    );
                                    di_j[i] = ClusterSequence::_bj_dij(&bj_jets[jet_i_idx], &bj_jets);
                                }
                                //if jet has NN of tail then used to be jetA
                                if jet_i_idx == end_idx {
                                    bj_jets[jet_i_idx].set_nn_jet(Some(min_idx));
                                }
                            }
                            None => {}
                        };
                    }
                }
            };
        }
    }

    
    fn do_jet_jet_recombination_step(&mut self, jet_a_idx: usize, jet_b_idx: usize, min_dij: f64) {
        // recombine jet_a and jet_b into new_jet
        let new_jet = self
            .jetdef
            .recombine(&self.particles[jet_a_idx], &self.particles[jet_b_idx]);

        // push new_jet onto particles
        self.particles.push(new_jet);

        // get it's index (new size of vec - 1)
        let new_idx = self.particles.len() - 1;

        // update history of jet_a and jet_b respectively
        let hist_a_idx = self.particles[jet_a_idx].cluster_hist_index();
        let hist_b_idx = self.particles[jet_b_idx].cluster_hist_index();

        self.add_recomb_to_history(hist_a_idx, hist_b_idx, new_idx, min_dij);
    }

    fn do_jet_beam_recombination_step(&mut self, jet_index: usize, d_ij: f64) {
        self.add_recomb_to_history(
            jet_index,
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

        self.particles[jet_index].set_cluster_hist_index(child_hist_index);
    }

    fn _bj_dij(jet: &BriefJet, jets: &[BriefJet]) -> f64 {
        let mut kt2 = jet.kt2();

        jet.nn_jet_index.map(|index| {
            if jets[index].kt2() < kt2 {
                kt2 = jets[index].kt2();
            }
        });

        jet.nn_dist * kt2
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
    fn bj_dist<J: ProxyJet>(jet_a: &J, jet_b: &J) -> f64 {
        let dphi: f64 = PI - f64::abs(PI - f64::abs(jet_a.phi() - jet_b.phi()));
        let deta: f64 = jet_a.eta() - jet_b.eta();
        dphi * dphi + deta * deta
    }

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
        let mut nn_dist = self.jetdef.r * self.jetdef.r;
        let mut nn: Option<usize> = None;

        // basically splitting for loop into two parts
        // can prolly use this in for loop with a || stop condition
        if head_idx < curr_idx {
            jets[head_idx..curr_idx]
                .iter()
                .enumerate()
                .for_each(|(jet_b_idx, jet_b)| {
                    let dist = ClusterSequence::bj_dist(&jets[curr_idx], jet_b);
                    if dist < nn_dist {
                        nn_dist = dist;
                        nn = Some(jet_b_idx);
                    }
                });
        }

        if tail_idx > curr_idx {
            jets[curr_idx + 1..=tail_idx]
                .iter()
                .enumerate()
                .for_each(|(jet_b_idx, jet_b)| {
                    let dist = ClusterSequence::bj_dist(&jets[curr_idx], jet_b);
                    if dist < nn_dist {
                        nn_dist = dist;
                        nn = Some(jet_b_idx + curr_idx + 1);
                    }
                });
        }
    }
    #[inline]
    fn bj_set_nn_crosscheck<'b, 'a, J: ProxyJet>(&'b self, jets: &'a mut [J]) -> () {
        let mut nn_dist = self.jetdef.r * self.jetdef.r;
        let mut nn: Option<usize> = None; //ENTIRE scope
        let n = jets.len();
        for i in 0..n {
            let dist = ClusterSequence::bj_dist(&jets[n], &jets[i]);
            if dist < nn_dist {
                nn_dist = dist;
                nn = Some(i); //HERE NN has ref to jets
            }
            if dist < jets[i].nn_dist() {
                jets[i].set_nn_dist(dist); //change jets &mut
                jets[i].set_nn_jet(Some(n));
            }
        }
        jets[n].set_nn_dist(nn_dist);
        jets[n].set_nn_jet(nn);
    }
}
