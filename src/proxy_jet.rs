#[derive(Debug, Clone)]
pub struct BriefJet {
    eta: f64,
    phi: f64,
    kt2: f64,
    pub nn_dist: f64, // either nn_dist == option, or only read dist if nn_jet != None
    pub nn_jet: Option<Box<BriefJet>>, //TODO: investigate if box is best way to do this
    index: usize,     // either index == option, or only read index if nn_jet != None
}

// pub struct TiledJet {
//     eta: f64,
//     phi: f64,
//     kt2: f64,
//     nn_dist: Option<f64>,// either nn_dist == option, or only read dist if nn_jet != None
//     nn_jet: Option<PseudoJet>,
//     index: Option<usize>, // either index == option, or only read index if nn_jet != None
// } //include linked lists for NN dist

impl BriefJet {
    fn new(eta: f64, phi: f64, kt2: f64, _r2: f64) -> Self {
        BriefJet {
            eta,
            phi,
            kt2,
            nn_dist: _r2,
            nn_jet: None,
            index: 0,
        }
    }
}

//for this trait include all generics
// include NN information
pub trait ProxyJet // problem is that jets do not have 4 mom only rap phi kt2 call this proxyjet
{
    fn eta(&self) -> f64;

    fn phi(&self) -> f64;

    fn kt2(&self) -> f64;

    //TODO: how to call different types of jet types
    fn nn_jet(&self) -> Option<Box<BriefJet>>
    where
        Self: Sized;

    fn nn_dist(&self) -> f64;

    // #[inline]
    // fn bj_diJ(&self) -> f64 {
    //     let kt2 = self.kt2();
    //     match self.NN {}
    // }
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

    #[inline]
    fn nn_jet(&self) -> Option<Box<BriefJet>> {
        self.nn_jet.clone()
    }
}
