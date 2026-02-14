use std::ops;

#[allow(non_snake_case)]
pub struct PseudoJet {
    _px: f64,
    _py: f64,
    _pz: f64,
    _E: f64,
    _rap: Option<f64>,
    _phi: Option<f64>,
    _kt2: Option<f64>, // _cluster_hist_index, _user_index (no clue what these are for yet)
}

#[allow(non_snake_case)]
impl PseudoJet {
    //TODO: investigate if unwraps is safe and does not panic crash all code

    const PSEUDOJET_INVALID_PHI: f64 = -100.0;
    const PSEUDOJET_INVALID_RAP: f64 = -1e200;

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
    pub fn rap(&self) -> f64 {
        self._rap.unwrap()
    }

    #[inline]
    pub fn phi(&self) -> f64 {
        self._phi.unwrap()
    }

    #[inline]
    pub fn kt2(&self) -> f64 {
        self._kt2.unwrap()
    }

    #[inline]
    pub fn pt2(&self) -> f64 {
        self._kt2.unwrap()
    }

    #[inline]
    pub fn pt(&self) -> f64 {
        self._kt2.unwrap().sqrt()
    }

    /// returns the squared transverse mass = kt^2+m^2
    #[inline]
    pub fn perp2(&self) -> f64 {
        self._kt2.unwrap()
    }

    /// returns the transverse mass = sqrt(kt^2+m^2)

    #[inline]
    pub fn perp(&self) -> f64 {
        self._kt2.unwrap().sqrt()
    }

    /// returns the squared transverse mass = kt^2+m^2
    #[inline]
    pub fn m2(&self) -> f64 {
        (self._E + self._pz) * (self._E - self._pz) - self._kt2.unwrap()
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
        self._kt2.unwrap() + self._pz * self._pz
    }

    /// return the 3-vector modulus = sqrt(px^2+py^2+pz^2)
    #[inline]
    pub fn modp(&self) -> f64 {
        self.modp2().sqrt()
    }

    /// return the transverse energy
    #[inline]
    pub fn et(&self) -> f64 {
        if self._kt2 == Some(0.0) {
            0.0
        } else {
            self._E / (1.0 + self._pz * self._pz / self._kt2.unwrap()).sqrt()
        }
    }
    /// return the transverse energy squared
    #[inline]
    pub fn et2(&self) -> f64 {
        if self._kt2 == Some(0.0) {
            0.0
        } else {
            self._E * self._E / (1.0 + self._pz * self._pz / self._kt2.unwrap())
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
        self._kt2 = Some(self.px() * self.px() + self.py() * self.py());
        self._phi = Some(Self::PSEUDOJET_INVALID_PHI);
        self._rap = Some(Self::PSEUDOJET_INVALID_RAP);
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

// implement operators as traits

impl ops::Add<PseudoJet> for PseudoJet {
    type Output = PseudoJet;

    fn add(self, other: PseudoJet) -> PseudoJet {
        PseudoJet {
            _px: self._px + other._px,
            _py: self._py + other._py,
            _pz: self._pz + other._pz,
            _E: self._E + other._E,
            _rap: None,
            _phi: None,
            _kt2: None,
        }
    }
}

// Summing a jet into this jet
impl ops::AddAssign<PseudoJet> for PseudoJet {
    fn add_assign(&mut self, other: PseudoJet) -> () {
        self._px += other._px;
        self._py += other._py;
        self._pz += other._pz;
        self._E += other._E;

        // now we need to reinit our values for this jet
        self.finish_init();
    }
}

impl ops::Sub<PseudoJet> for PseudoJet {
    type Output = PseudoJet;

    fn sub(self, other: PseudoJet) -> PseudoJet {
        PseudoJet {
            _px: self._px - other._px,
            _py: self._py - other._py,
            _pz: self._pz - other._pz,
            _E: self._E - other._E,
            _kt2: Some(self._kt2.unwrap() - other._kt2.unwrap()),
            _phi: None,
            _rap: None,
        }
    }
}

// Summing a jet into this jet
impl ops::SubAssign<PseudoJet> for PseudoJet {
    fn sub_assign(&mut self, other: PseudoJet) -> () {
        self._px -= other._px;
        self._py -= other._py;
        self._pz -= other._pz;
        self._E -= other._E;

        // now we need to reinit our values for this jet
        self.finish_init();
    }
}

impl ops::Mul<f64> for PseudoJet {
    type Output = PseudoJet;

    fn mul(self, scalar: f64) -> PseudoJet {
        // TODO: check if rap and phi is valid first (only needed for thread safety)
        //create new jet class
        let new_jet = PseudoJet {
            _px: self._px * scalar,
            _py: self._py * scalar,
            _pz: self._pz * scalar,
            _E: self._E * scalar,
            _kt2: Some(self._kt2.unwrap() * scalar * scalar),
            _phi: None,
            _rap: None,
        };

        new_jet
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
        let unwrap_kt2 = self._kt2.unwrap();
        self._kt2 = Some(unwrap_kt2 * scalar * scalar);
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

#[cfg(test)]
mod tests {
    use super::PseudoJet;

    fn approx_eq(a: f64, b: f64) {
        let diff = (a - b).abs();
        assert!(diff < 1e-12, "left={a}, right={b}, diff={diff}");
    }

    fn make_jet(px: f64, py: f64, pz: f64, e: f64) -> PseudoJet {
        let kt2 = px * px + py * py;
        PseudoJet {
            _px: px,
            _py: py,
            _pz: pz,
            _E: e,
            _rap: Some(-1e200),
            _phi: Some(-100.0),
            _kt2: Some(kt2),
        }
    }

    #[test]
    fn accessors_return_expected_components() {
        let jet = make_jet(3.0, 4.0, 12.0, 15.0);
        assert_eq!(jet.E(), 15.0);
        assert_eq!(jet.e(), 15.0);
        assert_eq!(jet.px(), 3.0);
        assert_eq!(jet.py(), 4.0);
        assert_eq!(jet.pz(), 12.0);
        assert_eq!(jet.rap(), -1e200);
        assert_eq!(jet.phi(), -100.0);
        assert_eq!(jet.kt2(), 25.0);
        assert_eq!(jet.pt2(), 25.0);
        assert_eq!(jet.perp2(), 25.0);
    }

    #[test]
    fn transverse_quantities_are_consistent() {
        let jet = make_jet(3.0, 4.0, 12.0, 15.0);
        assert_eq!(jet.pt(), 5.0);
        assert_eq!(jet.perp(), 5.0);
        assert_eq!(jet.m2(), 56.0);
        approx_eq(jet.m(), 56.0f64.sqrt());
        assert_eq!(jet.mperp2(), 81.0);
        assert_eq!(jet.mperp(), 9.0);
        assert_eq!(jet.mt2(), 81.0);
        assert_eq!(jet.mt(), 9.0);
    }

    #[test]
    fn three_vector_and_angle_quantities_are_consistent() {
        let jet = make_jet(3.0, 4.0, 12.0, 15.0);
        assert_eq!(jet.modp2(), 169.0);
        assert_eq!(jet.modp(), 13.0);
        approx_eq(jet.cos_theta(), 12.0 / 13.0);
        approx_eq(jet.theta(), (12.0 / 13.0f64).acos());
    }

    #[test]
    fn transverse_energy_is_correct_for_non_zero_kt2() {
        let jet = make_jet(3.0, 4.0, 12.0, 15.0);
        assert_eq!(jet.et(), 75.0 / 13.0);
        assert_eq!(jet.et2(), 5625.0 / 169.0);
    }

    #[test]
    fn transverse_energy_is_zero_when_kt2_is_zero() {
        let jet = make_jet(0.0, 0.0, 2.0, 5.0);
        assert_eq!(jet.et(), 0.0);
        assert_eq!(jet.et2(), 0.0);
    }

    #[test]
    fn index_returns_components() {
        let jet = make_jet(1.0, 2.0, 3.0, 4.0);
        assert_eq!(jet[0], 1.0);
        assert_eq!(jet[1], 2.0);
        assert_eq!(jet[2], 3.0);
        assert_eq!(jet[3], 4.0);
    }

    #[test]
    #[should_panic(expected = "Index out of bounds")]
    fn index_panics_out_of_bounds() {
        let jet = make_jet(1.0, 2.0, 3.0, 4.0);
        let _ = jet[4];
    }

    #[test]
    fn add_returns_componentwise_sum() {
        let a = make_jet(1.0, 2.0, 3.0, 4.0);
        let b = make_jet(0.5, 1.5, 2.5, 3.5);
        let c = a + b;
        assert_eq!(c._px, 1.5);
        assert_eq!(c._py, 3.5);
        assert_eq!(c._pz, 5.5);
        assert_eq!(c._E, 7.5);
        assert_eq!(c._kt2, None);
        assert_eq!(c._phi, None);
        assert_eq!(c._rap, None);
    }

    #[test]
    fn add_assign_updates_components_and_reinitializes_cached_values() {
        let mut a = make_jet(1.0, 2.0, 3.0, 4.0);
        let b = make_jet(0.5, 1.5, 2.5, 3.5);
        a += b;
        assert_eq!(a._px, 1.5);
        assert_eq!(a._py, 3.5);
        assert_eq!(a._pz, 5.5);
        assert_eq!(a._E, 7.5);
        assert_eq!(a._kt2, Some(14.5));
        assert_eq!(a._phi, Some(-100.0));
        assert_eq!(a._rap, Some(-1e200));
    }

    #[test]
    fn sub_returns_componentwise_difference() {
        let a = make_jet(1.0, 2.0, 3.0, 4.0);
        let b = make_jet(0.5, 1.5, 2.5, 3.5);
        let c = a - b;
        assert_eq!(c._px, 0.5);
        assert_eq!(c._py, 0.5);
        assert_eq!(c._pz, 0.5);
        assert_eq!(c._E, 0.5);
        assert_eq!(c._kt2, Some(2.5));
        assert_eq!(c._phi, None);
        assert_eq!(c._rap, None);
    }

    #[test]
    fn sub_assign_updates_components_and_reinitializes_cached_values() {
        let mut a = make_jet(1.0, 2.0, 3.0, 4.0);
        let b = make_jet(0.5, 1.5, 2.5, 3.5);
        a -= b;
        assert_eq!(a._px, 0.5);
        assert_eq!(a._py, 0.5);
        assert_eq!(a._pz, 0.5);
        assert_eq!(a._E, 0.5);
        assert_eq!(a._kt2, Some(0.5));
        assert_eq!(a._phi, Some(-100.0));
        assert_eq!(a._rap, Some(-1e200));
    }

    #[test]
    fn mul_scales_four_momentum_and_cached_kt2() {
        let jet = make_jet(1.0, 2.0, 3.0, 4.0);
        let out = jet * 3.0;
        assert_eq!(out._px, 3.0);
        assert_eq!(out._py, 6.0);
        assert_eq!(out._pz, 9.0);
        assert_eq!(out._E, 12.0);
        assert_eq!(out._kt2, Some(45.0));
        assert_eq!(out._phi, None);
        assert_eq!(out._rap, None);
    }

    #[test]
    fn mul_by_scalar_is_supported_from_left_side() {
        let jet = make_jet(1.0, 2.0, 3.0, 4.0);
        let out = 2.0 * jet;
        assert_eq!(out._px, 2.0);
        assert_eq!(out._py, 4.0);
        assert_eq!(out._pz, 6.0);
        assert_eq!(out._E, 8.0);
        assert_eq!(out._kt2, Some(20.0));
    }

    #[test]
    fn mul_assign_scales_in_place() {
        let mut jet = make_jet(1.0, 2.0, 3.0, 4.0);
        jet *= 2.0;
        assert_eq!(jet._px, 2.0);
        assert_eq!(jet._py, 4.0);
        assert_eq!(jet._pz, 6.0);
        assert_eq!(jet._E, 8.0);
        assert_eq!(jet._kt2, Some(20.0));
    }

    #[test]
    fn div_scales_by_inverse() {
        let jet = make_jet(2.0, 4.0, 6.0, 8.0);
        let out = jet / 2.0;
        assert_eq!(out._px, 1.0);
        assert_eq!(out._py, 2.0);
        assert_eq!(out._pz, 3.0);
        assert_eq!(out._E, 4.0);
        assert_eq!(out._kt2, Some(5.0));
    }

    #[test]
    fn div_assign_scales_by_inverse_in_place() {
        let mut jet = make_jet(2.0, 4.0, 6.0, 8.0);
        jet /= 2.0;
        assert_eq!(jet._px, 1.0);
        assert_eq!(jet._py, 2.0);
        assert_eq!(jet._pz, 3.0);
        assert_eq!(jet._E, 4.0);
        assert_eq!(jet._kt2, Some(5.0));
    }

    #[test]
    fn partial_eq_for_pseudojet_compares_four_momentum() {
        let a = make_jet(1.0, 2.0, 3.0, 4.0);
        let b = make_jet(1.0, 2.0, 3.0, 4.0);
        let c = make_jet(1.0, 2.0, 3.0, 5.0);
        assert!(a == b);
        assert!(a != c);
    }

    #[test]
    fn partial_eq_for_f64_supports_zero_only() {
        let zero = make_jet(0.0, 0.0, 0.0, 0.0);
        let non_zero = make_jet(1.0, 0.0, 0.0, 0.0);
        assert!(zero == 0.0);
        assert!(non_zero != 0.0);
    }

    #[test]
    #[should_panic(
        expected = "Comparing a PseudoJet with a non-zero constant (double) is not allowed"
    )]
    fn partial_eq_for_f64_panics_for_non_zero_constant() {
        let jet = make_jet(1.0, 2.0, 3.0, 4.0);
        let _ = jet == 1.0;
    }
}
