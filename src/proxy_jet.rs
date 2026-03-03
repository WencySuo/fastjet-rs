use crate::constants::PI;
use crate::pseudo_jet::PseudoJet;
use std::cell::RefCell;
use std::rc::Rc;

#[derive(Debug, Clone, PartialEq)]
pub struct BriefJet {
    pub eta: f64,
    pub phi: f64,
    pub kt2: f64,
    pub nn_dist: f64, // either nn_dist == option, or only read dist if nn_jet != None
    pub nn_jet_index: Option<usize>, //TODO: investigate if box is best way to do this
    pub particle_index: usize, // either index == option, or only read index if nn_jet != None
}

pub struct EEBriefJet {
    pub kt2: f64,
    pub nn_dist: f64, // either nn_dist == option, or only read dist if nn_jet != None
    pub nn_jet_index: Option<usize>,
    pub particle_index: usize, // either index == option, or only read index if nn_jet != None
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
#[derive(Default)]
pub struct TiledJet {
    pub eta: f64,
    pub phi: f64,
    pub kt2: f64,
    pub nn_dist: f64, // either nn_dist == option, or only read dist if nn_jet != None
    pub nn_jet_index: Option<usize>, //TODO: investigate if box is best way to do this
    pub particle_index: usize, // either index == option, or only read index if nn_jet != None
    pub bj_jet_index: usize,
    pub tile_index: usize,
    //Pointers to other jets
    pub prev_jet: Option<Rc<RefCell<TiledJet>>>,
    pub next_jet: Option<Rc<RefCell<TiledJet>>>,
}

#[derive(Clone, Default)]
pub struct Tile {
    pub begin_tiles: [usize; 9],
    pub surrounding_tiles: std::ops::Range<usize>,
    pub rh_tiles: std::ops::Range<usize>,
    pub head: Option<Rc<RefCell<TiledJet>>>, // need rc to make linked lists of jets
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

// impl TiledJet {
//     pub fn eta(&self) -> f64 {
//         self.eta
//     }

//     #[inline]
//     pub fn phi(&self) -> f64 {
//         self.phi
//     }
// }

impl ProxyJet for TiledJet {
    #[inline]
    fn eta(&self) -> f64 {
        self.eta
    }

    #[inline]
    fn phi(&self) -> f64 {
        self.phi
    }
    #[inline]
    fn kt2(&self) -> f64 {
        self.kt2
    }

    #[inline]
    fn nn_dist(&self) -> f64 {
        self.nn_dist
    }

    fn create_jet_type(&'_ self) -> JetType<'_> {
        JetType::TiledJetType(self)
    }

    #[inline]
    fn particle_index(&self) -> usize {
        self.particle_index
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
    fn _bj_dist(jet_a: &TiledJet, jet_b: &TiledJet) -> f64 {
        let dphi: f64 = PI - f64::abs(PI - f64::abs(jet_a.phi() - jet_b.phi()));
        let deta: f64 = jet_a.eta() - jet_b.eta();
        dphi * dphi + deta * deta
    }

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

        // use a min which defaults to max if not found maye bfaster?
        // jet.nn_dist()
        //     * jet.kt2().min(
        //         jet.nn_jet_index()
        //             .map(|index| jets[index].kt2())
        //             .unwrap_or(f64::MAX),
        //     )
    }

    fn _bj_set_jetinfo(index: usize, pseudo_jet: &PseudoJet, r2: f64, kt2: f64) -> TiledJet {
        TiledJet {
            eta: *pseudo_jet.rap(),
            phi: *pseudo_jet.phi(),
            kt2,
            nn_dist: r2,
            nn_jet_index: None,
            // at init paritlc_index and bj_jet_index should be the same
            particle_index: index,
            bj_jet_index: index,
            tile_index: 0, // default value
            prev_jet: None,
            next_jet: None,
        }
    }
}

//TODO remove the need for this class but needed for _bj_dij
impl ProxyJet for Rc<RefCell<TiledJet>> {
    #[inline]
    fn eta(&self) -> f64 {
        self.borrow().eta
    }

    #[inline]
    fn phi(&self) -> f64 {
        self.borrow().phi
    }

    #[inline]
    fn kt2(&self) -> f64 {
        self.borrow().kt2
    }

    #[inline]
    fn particle_index(&self) -> usize {
        self.borrow().particle_index
    }

    #[inline]
    fn nn_jet_index(&self) -> Option<usize> {
        self.borrow().nn_jet_index
    }

    #[inline]
    fn nn_dist(&self) -> f64 {
        self.borrow().nn_dist
    }

    #[inline]
    fn set_nn_dist(&mut self, dist: f64) {
        self.borrow_mut().nn_dist = dist;
    }

    #[inline]
    fn set_nn_jet(&mut self, jet_index: Option<usize>) {
        self.borrow_mut().nn_jet_index = jet_index;
    }

    #[inline]
    fn create_jet_type(&'_ self) -> JetType<'_> {
        // JetType wants a &'a TiledJet, so borrow then return ref.
        // IMPORTANT: This can't return a reference that outlives the borrow guard,
        // so if you truly need JetType here, you’ll likely want to refactor JetType
        // to own/carry Rc instead. If create_jet_type isn't used on tiled jets yet,
        // you can `unimplemented!()` for now.
        unimplemented!(
            "JetType::TiledJetType requires a &'a TiledJet; refactor JetType or avoid this for Rc<RefCell<_>>"
        );
    }

    #[inline]
    fn _bj_dist(jet_a: &Self, jet_b: &Self) -> f64 {
        let a = jet_a.borrow();
        let b = jet_b.borrow();
        TiledJet::_bj_dist(&a, &b)
    }

    fn _bj_set_jetinfo(index: usize, pseudo_jet: &PseudoJet, r2: f64, kt2: f64) -> Self {
        Rc::new(RefCell::new(TiledJet::_bj_set_jetinfo(
            index, pseudo_jet, r2, kt2,
        )))
    }

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
}

impl BriefJet {
    // TODO: removed since not used yet
    // fn new(eta: f64, phi: f64, kt2: f64, _r2: f64) -> Self {
    //     BriefJet {
    //         eta,
    //
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
    fn eta(&self) -> f64;

    fn phi(&self) -> f64;

    // problem is that jets do not have 4 mom only rap phi kt2 call this proxyjet
    fn kt2(&self) -> f64;

    fn particle_index(&self) -> usize;

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
        jet.nn_dist()
            * jet.kt2().min(
                jet.nn_jet_index()
                    .map(|index| jets[index].kt2())
                    .unwrap_or(f64::MAX),
            )
    }
    fn _bj_set_jetinfo(index: usize, pseudo_jet: &PseudoJet, r2: f64, kt2: f64) -> Self;
}

pub enum JetType<'a> {
    BriefJetType(&'a BriefJet),
    TiledJetType(&'a TiledJet),
    EEBriefJetType(&'a EEBriefJet),
}

// impl BriefJet {
//     #[inline]
//     fn eta(&self) -> f64 {
//         self.eta
//     }

//     #[inline]
//     fn phi(&self) -> f64 {
//         self.phi
//     }
// }

impl ProxyJet for BriefJet {
    #[inline]
    fn eta(&self) -> f64 {
        self.eta
    }

    #[inline]
    fn phi(&self) -> f64 {
        self.phi
    }

    #[inline]
    fn kt2(&self) -> f64 {
        self.kt2
    }

    #[inline]
    fn particle_index(&self) -> usize {
        self.particle_index
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
    fn _bj_dist(jet_a: &BriefJet, jet_b: &BriefJet) -> f64 {
        let dphi: f64 = PI - f64::abs(PI - f64::abs(jet_a.phi() - jet_b.phi()));
        let deta: f64 = jet_a.eta() - jet_b.eta();
        //dphi * dphi + deta * deta
        // modest 3% perfomance improvement versus prev impl
        dphi.mul_add(dphi, deta * deta)
    }

    fn _bj_set_jetinfo(index: usize, pseudo_jet: &PseudoJet, r2: f64, kt2: f64) -> BriefJet {
        BriefJet {
            eta: *pseudo_jet.rap(),
            phi: *pseudo_jet.phi(),
            kt2,
            particle_index: index,
            nn_dist: r2,
            nn_jet_index: None,
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
    pub fn _bj_dist_accurate(jet_a: &EEBriefJet, jet_b: &EEBriefJet) -> f64 {
        let mut dist =
            1.0 - jet_a.nx() * jet_b.nx() - jet_a.ny() * jet_b.ny() - jet_a.nz() * jet_b.nz();

        if dist * dist < f64::EPSILON {
            dist = (jet_a.ny() * jet_b.nz() - jet_b.ny() * jet_a.nz()).powf(2.0)
                + (jet_a.nz() * jet_b.nz() - jet_b.nz() * jet_a.nx()).powf(2.0)
                + (jet_a.nx() * jet_b.ny() - jet_b.nx() * jet_a.ny()).powf(2.0);
            return dist;
        }
        dist * 2.0
    }
}

impl ProxyJet for EEBriefJet {
    #[inline]
    fn eta(&self) -> f64 {
        panic!(" Should not be calling Phi in EE BJ")
    }

    #[inline]
    fn phi(&self) -> f64 {
        panic!(" Should not be calling Phi in EE BJ")
    }

    #[inline]
    fn kt2(&self) -> f64 {
        self.kt2
    }

    #[inline]
    fn particle_index(&self) -> usize {
        self.particle_index
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

    fn _bj_dist(jet_a: &EEBriefJet, jet_b: &EEBriefJet) -> f64 {
        let dist =
            1.0 - jet_a.nx() * jet_b.nx() - jet_a.ny() * jet_b.ny() - jet_a.nz() * jet_b.nz();
        dist * 2.0
    }

    fn _bj_set_jetinfo(index: usize, pseudo_jet: &PseudoJet, r2: f64, kt2: f64) -> EEBriefJet {
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
        EEBriefJet {
            kt2,
            particle_index: index,
            nn_dist: r2,
            nn_jet_index: None,
            nx,
            ny,
            nz,
        }
    }
}
