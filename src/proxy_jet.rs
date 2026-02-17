#[derive(Debug, Clone, PartialEq)]
pub struct BriefJet {
    pub eta: f64,
    pub phi: f64,
    pub kt2: f64,
    pub nn_dist: f64, // either nn_dist == option, or only read dist if nn_jet != None
    pub nn_jet_index: Option<usize>, //TODO: investigate if box is best way to do this
    pub _jets_index: usize,     // either index == option, or only read index if nn_jet != None
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
pub struct TiledJet<> {}

impl BriefJet {
    fn new(eta: f64, phi: f64, kt2: f64, _r2: f64) -> Self {
        BriefJet {
            eta,
            phi,
            kt2,
            nn_dist: _r2,
            nn_jet_index: None,
            _jets_index: 0,
        }
    }

    // fn get_jet_type(&self, jet_type: JetType) -> &BriefJet {
    //     match jet_type {
    //         JetType::BriefJetType(_) => self,
    //         _ => panic!("Unsupported jet type"),
    //     }
    // }

    //TODO: check if mut is necessary
    fn nn_jet<'a>(&self, jets: &'a mut[BriefJet]) -> Option<&'a mut BriefJet> {
        self.nn_jet_index.map(|index| &mut jets[index])
    }
}

//for this trait include all generics
// include NN information
pub trait ProxyJet {
    // problem is that jets do not have 4 mom only rap phi kt2 call this proxyjet
    fn eta(&self) -> f64;

    fn phi(&self) -> f64;

    fn kt2(&self) -> f64;

    fn eta_jet_type(&self, jet: & JetType) -> f64 {
        match jet {
            JetType::BriefJetType(jet) => jet.eta(),
            _ => panic!("Unsupported jet type"),
        }
    }

    fn phi_jet_type(&self, jet: & JetType) -> f64 {
        match jet {
            JetType::BriefJetType(jet) => jet.phi(),
            _ => panic!("Unsupported jet type"),
        }
    }

    fn kt2_jet_type(&self, jet: & JetType) -> f64 {
        match jet {
            JetType::BriefJetType(jet) => jet.kt2(),
            _ => panic!("Unsupported jet type"),
        }
    }

    fn create_jet_type<'a>(&'a self) -> JetType;


    //TODO: how to call different types of jet types
    fn nn_jet_index<'a>(&'a self) -> Option<usize>;

    fn nn_dist(&self) -> f64;

    fn set_nn_dist(&mut self, dist: f64);

    fn set_nn_jet(&mut self, jet_index: Option<usize>);

    // #[inline]
    // fn bj_diJ(&self) -> f64 {
    //     let kt2 = self.kt2();
    //     match self.NN {}
    // }
    //
}

enum JetType<'a> {
    BriefJetType(&'a BriefJet),
    TiledJetType( TiledJet),
}


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
    fn nn_dist(&self) -> f64 {
        self.nn_dist
    }

    fn create_jet_type(&self) -> JetType{
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
}
