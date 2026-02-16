use crate::constants::PI;
use crate::proxy_jet::{BriefJet, ProxyJet};
use crate::pseudo_jet::PseudoJet;
use std::ptr::eq;

// TODO: add automatic strategies for clustering like BEST: N2Plain, N2Tiling
pub struct ClusterSequence {
    pub particles: Vec<PseudoJet>,
    pub jetdef: JetDefinition,
}

pub struct JetDefinition {
    pub algorithm: Algorithm,
    pub r: f64,
}

// TODO: implement other algorithms
pub enum Algorithm {
    AntiKt,
    Cambridge,
}

impl ClusterSequence {
    //TODO: verify if this is right place for constant
    pub fn new(particles: Vec<PseudoJet>, jetdef: JetDefinition) -> Self {
        ClusterSequence { particles, jetdef }
    }

    // main driver loop
    pub fn simple_N2_cluster() {}

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
    fn bj_dij<J: ProxyJet>(jet: &J) -> f64 {
        let kt2 = jet.kt2();
        jet.nn_dist()
            * (match jet.nn_jet() {
                Some(jet_b) => {
                    if jet_b.kt2() < kt2 {
                        jet_b.kt2()
                    } else {
                        kt2
                    }
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
    fn bj_set_nn_crosscheck<J: ProxyJet>(
        &self,
        jets:&[BriefJet],
    ) -> () {
        
        // let mut nn_dist = self.jetdef.r * self.jetdef.r;
        // let mut jet_b = head;
        
        // for 

        // while !std::ptr::eq(jet_b, tail) {
        //     let dist = ClusterSequence::bj_dist(jet, jet_b);
        //     if dist < nn_dist {
        //         nn_dist = dist;
        //         jet.nn_jet = Some(Box::new(jet_b.clone()));
        //     }
        //     if dist < jet_b.nn_dist() {
        //         jet_b.nn_dist = dist;
        //         jet_b.nn_jet = Some(Box::new(jet_b.clone()));
        //     }
        // }
        // jet.nn_dist = nn_dist;
    }
}
