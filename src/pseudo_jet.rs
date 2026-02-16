use std::ops;

#[allow(non_snake_case)]
pub struct PseudoJet {
    _px: f64,
    _py: f64,
    _pz: f64,
    _E: f64,
    _rap: f64,
    _phi: f64,
    _kt2: f64, // _cluster_hist_index, _user_index (no clue what these are for yet)
}



// briefjet seems to be just the minimum amount of values from PseudoJet (aka without 3 mom) with NN information 

#[allow(non_snake_case)]
impl PseudoJet {
    pub const PSEUDOJET_INVALID_PHI: f64 = -100.0;
    pub const PSEUDOJET_INVALID_RAP: f64 = -1e200;

    pub fn new(px: f64, py: f64, pz: f64, E: f64) -> Self {
        PseudoJet {
            _px: px,
            _py: py,
            _pz: pz,
            _E: E,
            _kt2: (px * px + py * py),
            _phi: Self::PSEUDOJET_INVALID_PHI,
            _rap: Self::PSEUDOJET_INVALID_RAP,
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
    pub fn rap(& self) -> f64 {
        self._rap
    }

    #[inline]
    pub fn phi(& self) -> f64 {
        self._phi
    }

    #[inline]
    pub fn kt2(&self) -> f64 {
        self._kt2
    }

    #[inline]
    pub fn pt2(&self) -> f64 {
        self.kt2()
    }

    #[inline]
    pub fn pt(&self) -> f64 {
        self.kt2().sqrt()
    }

    /// returns the squared transverse mass = kt^2+m^2
    #[inline]
    pub fn perp2(&self) -> f64 {
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
            _ => self._E / (1.0 + self._pz * self._pz / self.kt2()).sqrt()
        }
    }
    /// return the transverse energy squared
    #[inline]
    pub fn et2(&self) -> f64 {
        match self.kt2() {
            0.0 => 0.0,
            _ => self._E * self._E / (1.0 + self._pz * self._pz / self.kt2())
        }
    }

    /// cos of the polar angle
    /// should we have: min(1.0,max(-1.0,_pz/sqrt(modp2())));
    #[inline]
    pub fn cos_theta(&self) -> f64 {
        return f64::min(1.0, f64::max(-1.0, self._pz / self.modp()));
    }

    /// polar angle
    #[inline]
    pub fn theta(&self) -> f64 {
        return self.cos_theta().acos();
    }

    fn finish_init(&mut self) {
        self._kt2 = self.px() * self.px() + self.py() * self.py();
        self._phi = Self::PSEUDOJET_INVALID_PHI;
        self._rap = Self::PSEUDOJET_INVALID_RAP;
    }
}

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
fn add_assign(&mut self, other: PseudoJet) -> () {
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
            self._E - other._E
        )
    }
}

// Summing a jet into this jet
impl ops::SubAssign<PseudoJet> for PseudoJet {
    fn sub_assign(&mut self, other: PseudoJet) -> () {
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
    fn mul_assign(&mut self, scalar: f64) -> () {
        self._px *= scalar;
        self._py *= scalar;
        self._pz *= scalar;
        self._E *= scalar;
        self._kt2 = self.kt2() * scalar * scalar;
    }
}

impl ops::Div<f64> for PseudoJet {
    type Output = PseudoJet;

    fn div(self, rhs: f64) -> PseudoJet {
        self * (1.0 / rhs)
    }
}

impl ops::DivAssign<f64> for PseudoJet {
    fn div_assign(&mut self, scalar: f64) -> () {
        *self *= 1.0 / scalar;
    }
}

impl PartialEq<PseudoJet> for PseudoJet {
    fn eq(&self, other: &PseudoJet) -> bool {
        if self._px != other._px
            || self._py != other._py
            || self._pz != other._pz
            || self._E != other._E
        {
            return false;
        } else {
            return true;
        }
    }
}

impl PartialEq<f64> for PseudoJet {
    fn eq(&self, val: &f64) -> bool {
        match val {
            0.0 => return self._px == 0.0 && self._py == 0.0 && self._pz == 0.0 && self._E == 0.0,
            _ => panic!("Comparing a PseudoJet with a non-zero constant (double) is not allowed"),
        }
    }
}
