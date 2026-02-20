use crate::cluster_sequence::INVALID;
use crate::constants::PI;
use std::cell::OnceCell;
use std::ops;

#[allow(non_snake_case)]
#[derive(Clone)]
pub struct PseudoJet {
    _px: f64,
    _py: f64,
    _pz: f64,
    _E: f64,
    _rap: OnceCell<f64>,
    _phi: OnceCell<f64>,
    _kt2: OnceCell<f64>, // _cluster_hist_index, _user_index (no clue what these are for yet)
    _cluster_hist_index: usize,
}

// briefjet seems to be just the minimum amount of values from PseudoJet (aka without 3 mom) with NN information

#[allow(non_snake_case)]
impl PseudoJet {
    pub fn new(px: f64, py: f64, pz: f64, E: f64) -> Self {
        PseudoJet {
            _px: px,
            _py: py,
            _pz: pz,
            _E: E,
            _kt2: OnceCell::from(px * px + py * py),
            _phi: OnceCell::new(),
            _rap: OnceCell::new(),
            _cluster_hist_index: INVALID,
        }
    }

    #[inline]
    pub fn E(&self) -> f64 {
        self._E
    }

    #[inline]
    pub fn e(&self) -> f64 {
        self._E
    }

    #[inline]
    pub fn px(&self) -> f64 {
        self._px
    }

    #[inline]
    pub fn py(&self) -> f64 {
        self._py
    }

    #[inline]
    pub fn pz(&self) -> f64 {
        self._pz
    }

    #[inline]
    pub fn rap(&self) -> &f64 {
        self._rap.get_or_init(|| {
            let max_rap: f64 = 1e5;
            let mut rap;
            let kt2 = *self.kt2();
            if self.e() == self.pz().abs() && kt2 == 0.0 {
                let max_rap_here = max_rap + self.pz().abs();
                if self.pz() > 0.0 {
                    rap = max_rap_here;
                } else {
                    rap = -max_rap_here;
                }
            } else {
                let effective_m2 = self.m2().max(0.0);
                let e_plus_pz = self.e() + self.pz().abs();
                rap = 0.5 * ((kt2 + effective_m2) / (e_plus_pz * e_plus_pz)).ln();
                if self.pz() > 0.0 {
                    rap = -rap;
                }
            }
            rap
        })
    }

    #[inline]
    pub fn phi(&self) -> &f64 {
        self._phi.get_or_init(|| {
            let mut phi;
            if *(self.kt2()) == 0.0 {
                phi = 0.0;
            } else {
                phi = self.py().atan2(self.px())
            }

            if phi < 0.0 {
                phi += 2.0 * PI;
            }

            if phi >= 2.0 * PI {
                phi -= 2.0 * PI;
            }
            phi
        })
    }

    // #[inline]
    // pub fn _set_rap_phi() {
    //     let max_rap: f64 = 1e5;
    //     if self.kt2() == 0.0 {
    //         self._phi = 0.0;
    //     } else {
    //         self._phi = self.py().atan2(self.px())
    //     }

    //     if self._phi < 0.0 {
    //         self._phi += 2.0 * PI;
    //     }

    //     if self._phi >= 2.0 * PI {
    //         self._phi -= 2.0 * PI;
    //     }

    //     if self.e() == self.pz().abs() && self.kt2() == 0.0 {
    //         let max_rap_here = max_rap + self.pz().abs();
    //         if self.pz() > 0.0 {
    //             self._rap = max_rap_here;
    //         } else {
    //             self._rap = -max_rap_here;
    //         }
    //     } else {
    //         let effective_m2 = self.m2().max(0.0);
    //         let e_plus_pz = self.e() + self.pz().abs();
    //         self._rap = 0.5 * ((self.kt2() + effective_m2) / (e_plus_pz * e_plus_pz)).ln();
    //         if self.pz() > 0.0 {
    //             self._rap = -self._rap;
    //         }
    //     }
    // }

    #[inline]
    pub fn kt2(&self) -> &f64 {
        self._kt2
            .get_or_init(|| self._px * self._px + self._py * self._py)
    }

    #[inline]
    pub fn pt2(&self) -> &f64 {
        self.kt2()
    }

    #[inline]
    pub fn pt(&self) -> f64 {
        self.kt2().sqrt()
    }

    /// returns the squared transverse mass = kt^2+m^2
    #[inline]
    pub fn perp2(&self) -> &f64 {
        self.kt2()
    }

    /// returns the transverse mass = sqrt(kt^2+m^2)
    #[inline]
    pub fn perp(&self) -> f64 {
        self.kt2().sqrt()
    }

    /// returns the squared transverse mass = kt^2+m^2
    #[inline]
    pub fn m2(&self) -> f64 {
        (self._E + self._pz) * (self._E - self._pz) - self.kt2()
    }

    /// returns the transverse mass = sqrt(kt^2+m^2)
    #[inline]
    pub fn m(&self) -> f64 {
        // taken literally from CLHEP
        let mm: f64 = self.m2();
        if mm < 0.0 { -((-mm).sqrt()) } else { mm.sqrt() }
    }

    /// returns the squared transverse mass = kt^2+m^2
    #[inline]
    pub fn mperp2(&self) -> f64 {
        (self._E + self._pz) * (self._E - self._pz)
    }

    /// returns the transverse mass = sqrt(kt^2+m^2)
    #[inline]
    pub fn mperp(&self) -> f64 {
        self.mperp2().abs().sqrt()
    }

    /// returns the squared transverse mass = kt^2+m^2
    #[inline]
    pub fn mt2(&self) -> f64 {
        (self._E + self._pz) * (self._E - self._pz)
    }

    /// returns the transverse mass = sqrt(kt^2+m^2)
    #[inline]
    pub fn mt(&self) -> f64 {
        self.mperp2().abs().sqrt()
    }

    /// return the squared 3-vector modulus = px^2+py^2+pz^2
    #[inline]
    pub fn modp2(&self) -> f64 {
        self.kt2() + self._pz * self._pz
    }

    /// return the 3-vector modulus = sqrt(px^2+py^2+pz^2)
    #[inline]
    pub fn modp(&self) -> f64 {
        self.modp2().sqrt()
    }

    /// return the transverse energy
    #[inline]
    pub fn et(&self) -> f64 {
        match self.kt2() {
            0.0 => 0.0,
            _ => self._E / (1.0 + self._pz * self._pz / self.kt2()).sqrt(),
        }
    }
    /// return the transverse energy squared
    #[inline]
    pub fn et2(&self) -> f64 {
        match self.kt2() {
            0.0 => 0.0,
            _ => self._E * self._E / (1.0 + self._pz * self._pz / self.kt2()),
        }
    }

    /// cos of the polar angle
    /// should we have: min(1.0,max(-1.0,_pz/sqrt(modp2())));
    #[inline]
    pub fn cos_theta(&self) -> f64 {
        (self._pz / self.modp()).clamp(-1.0, 1.0)
    }

    /// polar angle
    #[inline]
    pub fn theta(&self) -> f64 {
        self.cos_theta().acos()
    }

    pub fn set_cluster_hist_index(&mut self, index: usize) {
        self._cluster_hist_index = index;
    }

    pub fn cluster_hist_index(&self) -> usize {
        self._cluster_hist_index
    }

    fn finish_init(&mut self) {
        self._kt2 = OnceCell::from(self.px() * self.px() + self.py() * self.py());
        self._phi = OnceCell::new();
        self._rap = OnceCell::new();
    }

    // TODO: investigate if we need custom sorting implementation for performance reasons
    pub fn sorted_by_pt<'a>(jets: &'a mut Vec<&'a PseudoJet>) -> &'a mut Vec<&'a PseudoJet> {
        jets.sort_by(|a, b| (-a.kt2()).total_cmp(&(-b.kt2())));
        jets
    }
}

// TODO: investigate if we ever need to compare raw PseudoJet objs as a whole and not just sort by fields
// impl Ord for PseudoJet {
//     fn cmp(&self, other: &Self) -> std::cmp::Ordering {
//         self._px.total_cmp(&other._px)
//             .then(self._py.total_cmp(&other._py))
//             .then(self._pz.total_cmp(&other._pz))
//             .then(self._E.total_cmp(&other._E))
//             .then(self._kt2.total_cmp(&other._kt2))
//             .then(self._rap.total_cmp(&other._rap))
//             .then(self._phi.total_cmp(&other._phi))
//             .then(self._cluster_hist_index.cmp(&other._cluster_hist_index))
//     }
// }

// impl PartialOrd for PseudoJet {
//     fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
//         Some(self.cmp(other))
//     }
// }

// impl Eq for PseudoJet {}

// impl PartialEq for PseudoJet {
//     fn eq(&self, other: &Self) -> bool {
//         self._px == other._px
//         && self._py == other._py
//         && self._pz == other._pz
//         && self._E == other._E
//         && self._rap == other._rap
//         && self._phi == other._phi
//         && self._kt2 == other._kt2
//         && self._cluster_hist_index == other._cluster_hist_index
//     }
// }

impl ops::Index<usize> for PseudoJet {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self._px,
            1 => &self._py,
            2 => &self._pz,
            3 => &self._E,
            _ => panic!("Index out of bounds"),
        }
    }
}

impl ops::Add<PseudoJet> for PseudoJet {
    type Output = PseudoJet;

    fn add(self, other: PseudoJet) -> PseudoJet {
        PseudoJet::new(
            self._px + other._px,
            self._py + other._py,
            self._pz + other._pz,
            self._E + other._E,
        )
    }
}

impl ops::Add<&PseudoJet> for &PseudoJet {
    type Output = PseudoJet;

    fn add(self, other: &PseudoJet) -> PseudoJet {
        PseudoJet::new(
            self._px + other._px,
            self._py + other._py,
            self._pz + other._pz,
            self._E + other._E,
        )
    }
}

// Summing a jet into this jet
impl ops::AddAssign<PseudoJet> for PseudoJet {
    fn add_assign(&mut self, other: PseudoJet) {
        self._px += other._px;
        self._py += other._py;
        self._pz += other._pz;
        self._E += other._E;

        self.finish_init();
    }
}

impl ops::Sub<PseudoJet> for PseudoJet {
    type Output = PseudoJet;

    fn sub(self, other: PseudoJet) -> PseudoJet {
        PseudoJet::new(
            self._px - other._px,
            self._py - other._py,
            self._pz - other._pz,
            self._E - other._E,
        )
    }
}

// Summing a jet into this jet
impl ops::SubAssign<PseudoJet> for PseudoJet {
    fn sub_assign(&mut self, other: PseudoJet) {
        self._px -= other._px;
        self._py -= other._py;
        self._pz -= other._pz;
        self._E -= other._E;

        self.finish_init();
    }
}

impl ops::Mul<f64> for PseudoJet {
    type Output = PseudoJet;

    fn mul(self, scalar: f64) -> PseudoJet {
        // TODO: check if rap and phi is valid first (only needed for thread safety)
        PseudoJet::new(
            self._px * scalar,
            self._py * scalar,
            self._pz * scalar,
            self._E * scalar,
        )
    }
}

impl ops::Mul<PseudoJet> for f64 {
    type Output = PseudoJet;

    fn mul(self, scalar: PseudoJet) -> PseudoJet {
        scalar * self
    }
}

impl ops::MulAssign<f64> for PseudoJet {
    fn mul_assign(&mut self, scalar: f64) {
        self._px *= scalar;
        self._py *= scalar;
        self._pz *= scalar;
        self._E *= scalar;
        self._kt2 = OnceCell::new();
        self._rap = OnceCell::new();
        self._phi = OnceCell::new();
    }
}

impl ops::Div<f64> for PseudoJet {
    type Output = PseudoJet;

    fn div(self, rhs: f64) -> PseudoJet {
        self * (1.0 / rhs)
    }
}

impl ops::DivAssign<f64> for PseudoJet {
    fn div_assign(&mut self, scalar: f64) {
        *self *= 1.0 / scalar;
    }
}

impl PartialEq<PseudoJet> for PseudoJet {
    fn eq(&self, other: &PseudoJet) -> bool {
        !(self._px != other._px
            || self._py != other._py
            || self._pz != other._pz
            || self._E != other._E)
    }
}

impl PartialEq<f64> for PseudoJet {
    fn eq(&self, val: &f64) -> bool {
        match val {
            0.0 => self._px == 0.0 && self._py == 0.0 && self._pz == 0.0 && self._E == 0.0,
            _ => panic!("Comparing a PseudoJet with a non-zero constant (double) is not allowed"),
        }
    }
}
