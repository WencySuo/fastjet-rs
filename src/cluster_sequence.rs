use crate::constants::PI;
use crate::proxy_jet::BriefJet;
use crate::proxy_jet::ProxyJet;
use crate::pseudo_jet::PseudoJet;
use std::mem; 

// TODO: add automatic strategies for clustering like BEST: N2Plain, N2Tiling
pub struct ClusterSequence {
    pub particles: Vec<PseudoJet>,
    pub jetdef: JetDefinition,
    _invr2: f64,
    r2: f64,
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
        JetDefinition { algorithm, r, scheme }
    }
    
    //TODO: should add reset function to make it faster
    pub fn recombine(&self, jet_a: &PseudoJet, jet_b: &PseudoJet) -> PseudoJet {
        match self.scheme {
            RecombinationScheme::EScheme => {
                let jet_ab = jet_b + jet_a; //keeps mem in jet_b
                jet_ab
            }, 
            _ => {
                panic!("scheme not yet supported")
            }
        }
    }
}

impl ClusterSequence {
    //TODO: verify if this is right place for constant
    pub fn new(particles: Vec<PseudoJet>, jetdef: JetDefinition) -> Self {
        let r2 = jetdef.r * jetdef.r;
        let _invr2 = 1.0 / r2;
        ClusterSequence { particles, jetdef, _invr2, r2 }
    }

    // main driver loop
    pub fn simple_n2_cluster(&self) {
        let n = self.particles.len();

        //
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
            *jet += self._bj_dij(&bj_jets[i], &bj_jets);
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
            let end_index = n - 1 - i;
            
            // impossible to not have min deltaR
            let (min_index, min_dij) = di_j[0..=end_index].iter().enumerate().min_by(|a, b| a.1.total_cmp(&b.1)).unwrap();
            
            let mut jet_a = &bj_jets[min_index];
            let mut jet_b: Option<&BriefJet> = jet_a.nn_jet_index.map(|index| &bj_jets[index]);
            
            // normalize
            let min_dij = min_dij * self._invr2;
            
            //check if jet_a has neighbor
            match jet_b {
                Some(jet) => {
                     // nn is the new index we want to store the recombined index at 
                     //modify jetB for BJ with new info
                    do_jet_jet_recombination_step(jet_a, jet_b, min_dij); 
                }
                None => {
                    // jet beam recombination 
                    do_jet_beam_recombination_step();
                }
            }
            
            //swap jetA with tail vector
            bj_jets.swap(min_index, end_index);
            // now swap the min distances loc also
            di_j.swap(min_index, end_index);
            
            // update the diJ[] 
            
            
            
        }
    }
    
    fn do_jet_jet_recombination_step(&self, jet_a_idx: &mut usize, jet_b_idx: &mut usize, min_dij: f64) {
        // call recombiner to recombine 
        self.jetdef.recombine(self.particles[*jet_a_idx], self.particles[*jet_b_idx]);
    }

    fn do_jet_beam_recombination_step() {
        // jet beam recombination 

    }

    fn _bj_dij(&self, jet: &BriefJet, jets: &[BriefJet]) -> f64 {
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
