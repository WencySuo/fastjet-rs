use crate::constants::PI;
use crate::pseudo_jet::PseudoJet;
use std::rc::Rc;

#[derive(Debug, Clone, PartialEq)]
pub struct BriefJet {
    pub eta: f64,
    pub phi: f64,
    pub kt2: f64,
    pub nn_dist: f64, // either nn_dist == option, or only read dist if nn_jet != None
    pub nn_jet_index: Option<usize>, //TODO: investigate if box is best way to do this
    pub _jets_index: usize, // either index == option, or only read index if nn_jet != None
}

pub struct EEBriefJet {
    pub kt2: f64,
    pub nn_dist: f64, // either nn_dist == option, or only read dist if nn_jet != None
    pub nn_jet_index: Option<usize>,
    pub _jets_index: usize, // either index == option, or only read index if nn_jet != None
    pub nx: f64,
    pub ny: f64,
    pub nz: f64,
}

// pub struct TiledJet {
//     eta: f64,
//     phi: f64,
//     kt2: f64,
//     nn_dist: Option<f64>,// either nn_dist == option, or only read dist if nn_jet != None
//     nn_jet: Option<PseudoJet>,
//     index: Option<usize>, // either index == option, or only read index if nn_jet != None
// } //include linked lists for NN dist

//empty struct for now
pub struct TiledJet {
    pub eta: f64,
    pub phi: f64,
    pub kt2: f64,
    pub nn_dist: f64, // either nn_dist == option, or only read dist if nn_jet != None
    pub nn_jet_index: Option<usize>, //TODO: investigate if box is best way to do this
    pub _jets_index: usize, // either index == option, or only read index if nn_jet != None
    pub _tile_index: usize,
    pub dij_posn: usize,
    //Pointers to other Tiles
    pub prev_tile: Rc<TiledJet>,
    pub next_tile: Rc<TiledJet>,
}

#[derive(Clone)]
pub struct Tile {
    pub begin_tiles: [usize; 9],
    pub surrounding_tiles: std::ops::Range<usize>,
    pub rh_tiles: std::ops::Range<usize>,
    pub head: Option<usize>, // index of TiledJet
}

impl Tile {
    pub fn new() -> Self {
        Tile {
            begin_tiles: [0; 9],
            surrounding_tiles: 0..0,
            rh_tiles: 0..0,
            head: None,
        }
    }
}

impl BriefJet {
    // TODO: removed since not used yet
    // fn new(eta: f64, phi: f64, kt2: f64, _r2: f64) -> Self {
    //     BriefJet {
    //         eta,
    //         phi,
    //         kt2,
    //         nn_dist: _r2,
    //         nn_jet_index: None,
    //         _jets_index: 0,
    //     }
    // }

    // fn get_jet_type(&self, jet_type: JetType) -> &BriefJet {
    //     match jet_type {
    //         JetType::BriefJetType(_) => self,
    //         _ => panic!("Unsupported jet type"),
    //     }
    // }

    // TODO: removed since not used yet
    // //TODO: check if mut is necessary
    // fn nn_jet<'a>(&self, jets: &'a mut[BriefJet]) -> Option<&'a mut BriefJet> {
    //     self.nn_jet_index.map(|index| &mut jets[index])
    // }
}

//for this trait include all generics
// include NN information
pub trait ProxyJet {
    // problem is that jets do not have 4 mom only rap phi kt2 call this proxyjet
    fn kt2(&self) -> f64;

    fn jets_index(&self) -> usize;

    fn eta_jet_type(&self, jet: &JetType) -> f64 {
        match jet {
            JetType::BriefJetType(jet) => jet.eta(),
            _ => panic!("Unsupported jet type"),
        }
    }

    fn phi_jet_type(&self, jet: &JetType) -> f64 {
        match jet {
            JetType::BriefJetType(jet) => jet.phi(),
            _ => panic!("Unsupported jet type"),
        }
    }

    fn kt2_jet_type(&self, jet: &JetType) -> f64 {
        match jet {
            JetType::BriefJetType(jet) => jet.kt2(),
            _ => panic!("Unsupported jet type"),
        }
    }

    fn create_jet_type(&'_ self) -> JetType<'_>;

    //TODO: how to call different types of jet types
    fn nn_jet_index(&self) -> Option<usize>;

    fn nn_dist(&self) -> f64;

    fn set_nn_dist(&mut self, dist: f64);

    fn set_nn_jet(&mut self, jet_index: Option<usize>);

    fn _bj_dist(jet_a: &Self, jet_b: &Self) -> f64;

    #[inline]
    fn _bj_dij<J: ProxyJet>(jet: &J, jets: &[J]) -> f64 {
        let mut kt2 = jet.kt2();
        if let Some(index) = jet.nn_jet_index() {
            let kt2_b = jets[index].kt2();
            if kt2_b < kt2 {
                kt2 = kt2_b;
            }
        }
        jet.nn_dist() * kt2
    }

    fn _bj_set_jetinfo(index:usize, pseudo_jet: &PseudoJet, r2:f64, kt2:f64) -> Self;
}

pub enum JetType<'a> {
    BriefJetType(&'a BriefJet),
    TiledJetType(TiledJet),
    EEBriefJetType(&'a EEBriefJet),
}

impl BriefJet {
    #[inline]
    fn eta(&self) -> f64 {
        self.eta
    }

    #[inline]
    fn phi(&self) -> f64 {
        self.phi
    }
}

impl ProxyJet for BriefJet {
    #[inline]
    fn kt2(&self) -> f64 {
        self.kt2
    }

    #[inline]
    fn jets_index(&self) -> usize {
        self._jets_index
    }

    #[inline]
    fn nn_dist(&self) -> f64 {
        self.nn_dist
    }

    fn create_jet_type(&'_ self) -> JetType<'_> {
        JetType::BriefJetType(self)
    }

    #[inline]
    fn nn_jet_index(&self) -> Option<usize> {
        self.nn_jet_index
    }

    #[inline]
    fn set_nn_dist(&mut self, dist: f64) {
        self.nn_dist = dist;
    }

    #[inline]
    fn set_nn_jet(&mut self, jet_index: Option<usize>) {
        self.nn_jet_index = jet_index;
    }
    
    #[inline]
    fn _bj_dist(jet_a:&BriefJet, jet_b: &BriefJet) -> f64 {
        let dphi: f64 = PI - f64::abs(PI - f64::abs(jet_a.phi() - jet_b.phi()));
        let deta: f64 = jet_a.eta() - jet_b.eta();
        dphi * dphi + deta * deta
    }

    fn _bj_set_jetinfo(index:usize, pseudo_jet: &PseudoJet, r2:f64, kt2:f64) -> BriefJet {
        BriefJet{
            eta: *pseudo_jet.rap(),
            phi: *pseudo_jet.phi(),
            kt2,
            _jets_index: index,
            nn_dist: r2,
            nn_jet_index: None
        }
    }

}

impl EEBriefJet {
    #[inline]
    fn nx(&self) -> f64 {
        self.nx
    }
    
    #[inline]
    fn ny(&self) -> f64 {
        self.ny
    }
    
    #[inline]
    fn nz(&self) -> f64 {
        self.nz
    }
    
    // method derived from PanScales for more accurate distance calculation
    pub fn _bj_dist_accurate(jet_a:&EEBriefJet, jet_b:&EEBriefJet) -> f64 {
        let mut dist = 1.0
            - jet_a.nx() * jet_b.nx()
            - jet_a.ny() * jet_b.ny()
            - jet_a.nz() * jet_b.nz();
        
        if dist*dist < f64::EPSILON {
            dist = (jet_a.ny() * jet_b.nz() - jet_b.ny() * jet_a.nz()).powf(2.0) 
                + (jet_a.nz() * jet_b.nz() - jet_b.nz() * jet_a.nx()).powf(2.0)
                + (jet_a.nx() * jet_b.ny() - jet_b.nx() * jet_a.ny()).powf(2.0);
            return dist; 
        }
        dist * 2.0
    }
}

impl ProxyJet for EEBriefJet{
    #[inline]
    fn kt2(&self) -> f64 {
        self.kt2
    }

    #[inline]
    fn jets_index(&self) -> usize {
        self._jets_index
    }
    
    fn create_jet_type(&'_ self) -> JetType<'_> {
        JetType::EEBriefJetType(self)
    }
    
    #[inline]
    fn nn_dist(&self) -> f64 {
        self.nn_dist
    }

    fn nn_jet_index(&self) -> Option<usize> {
        self.nn_jet_index
    }
    
    #[inline]
    fn set_nn_dist(&mut self, dist: f64) {
        self.nn_dist = dist;
    }

    #[inline]
    fn set_nn_jet(&mut self, jet_index: Option<usize>) {
        self.nn_jet_index = jet_index;
    }

    fn _bj_dist(jet_a:&EEBriefJet, jet_b:&EEBriefJet) -> f64 {
        let dist = 1.0
            - jet_a.nx() * jet_b.nx()
            - jet_a.ny() * jet_b.ny()
            - jet_a.nz() * jet_b.nz();
        dist * 2.0
    }

    fn _bj_set_jetinfo(index:usize, pseudo_jet: &PseudoJet, r2:f64, kt2:f64) -> EEBriefJet {
        let mut norm = pseudo_jet.modp2();
        let mut nx = 0.0;
        let mut ny = 0.0;
        let mut nz = 1.0;
        if norm > 0.0 {
            norm = 1.0 / norm.sqrt();
            nx = pseudo_jet.px() * norm;
            ny = pseudo_jet.py() * norm;
            nz = pseudo_jet.pz() * norm;
        }
        EEBriefJet{
            kt2,
            _jets_index: index,
            nn_dist: r2,
            nn_jet_index: None,
            nx,
            ny,
            nz
        }
    }
}
