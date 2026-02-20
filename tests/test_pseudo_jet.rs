#[cfg(test)]
mod tests {
    #[cfg(test)]
    mod tests {
        use fastjet_rs::pseudo_jet::PseudoJet;

        fn approx_eq(a: f64, b: f64) {
            let diff = (a - b).abs();
            assert!(diff < 1e-12, "left={a}, right={b}, diff={diff}");
        }

        #[test]
        fn accessors_return_expected_components() {
            let jet = PseudoJet::new(3.0, 4.0, 12.0, 15.0);
            assert_eq!(jet.E(), 15.0);
            assert_eq!(jet.e(), 15.0);
            assert_eq!(jet.px(), 3.0);
            assert_eq!(jet.py(), 4.0);
            assert_eq!(jet.pz(), 12.0);
            assert_eq!(*jet.rap(), 1.0986122886681098);
            assert_eq!(*jet.phi(), 0.9272952180016122);
            assert_eq!(*jet.kt2(), 25.0);
            assert_eq!(*jet.pt2(), 25.0);
            assert_eq!(*jet.perp2(), 25.0);
        }

        #[test]
        fn transverse_quantities_are_consistent() {
            let jet = PseudoJet::new(3.0, 4.0, 12.0, 15.0);
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
            let jet = PseudoJet::new(3.0, 4.0, 12.0, 15.0);
            assert_eq!(jet.modp2(), 169.0);
            assert_eq!(jet.modp(), 13.0);
            approx_eq(jet.cos_theta(), 12.0 / 13.0);
            approx_eq(jet.theta(), (12.0 / 13.0f64).acos());
        }

        #[test]
        fn transverse_energy_is_correct_for_non_zero_kt2() {
            let jet = PseudoJet::new(3.0, 4.0, 12.0, 15.0);
            assert_eq!(jet.et(), 75.0 / 13.0);
            assert_eq!(jet.et2(), 5625.0 / 169.0);
        }

        #[test]
        fn transverse_energy_is_zero_when_kt2_is_zero() {
            let jet = PseudoJet::new(0.0, 0.0, 2.0, 5.0);
            assert_eq!(jet.et(), 0.0);
            assert_eq!(jet.et2(), 0.0);
        }

        #[test]
        fn index_returns_components() {
            let jet = PseudoJet::new(1.0, 2.0, 3.0, 4.0);
            assert_eq!(jet[0], 1.0);
            assert_eq!(jet[1], 2.0);
            assert_eq!(jet[2], 3.0);
            assert_eq!(jet[3], 4.0);
        }

        #[test]
        #[should_panic(expected = "Index out of bounds")]
        fn index_panics_out_of_bounds() {
            let jet = PseudoJet::new(1.0, 2.0, 3.0, 4.0);
            let _ = jet[4];
        }

        #[test]
        fn add_returns_componentwise_sum() {
            let a = PseudoJet::new(1.0, 2.0, 3.0, 4.0);
            let b = PseudoJet::new(0.5, 1.5, 2.5, 3.5);
            let c = a + b;
            assert_eq!(c.px(), 1.5);
            assert_eq!(c.py(), 3.5);
            assert_eq!(c.pz(), 5.5);
            assert_eq!(c.e(), 7.5);
            assert_eq!(*c.kt2(), 14.5);
            assert_eq!(*c.phi(), 1.1659045405098132);
            assert_eq!(*c.rap(), 0.9359010884507957);
        }

        #[test]
        fn add_assign_updates_components_and_reinitializes_cached_values() {
            let mut a = PseudoJet::new(1.0, 2.0, 3.0, 4.0);
            let b = PseudoJet::new(0.5, 1.5, 2.5, 3.5);
            a += b;
            assert_eq!(a.px(), 1.5);
            assert_eq!(a.py(), 3.5);
            assert_eq!(a.pz(), 5.5);
            assert_eq!(a.e(), 7.5);
            assert_eq!(*a.kt2(), 14.5);
            assert_eq!(*a.phi(), 1.1659045405098132);
            assert_eq!(*a.rap(), 0.9359010884507957);
        }

        #[test]
        fn sub_returns_componentwise_difference() {
            let a = PseudoJet::new(1.0, 2.0, 3.0, 4.0);
            let b = PseudoJet::new(0.5, 1.5, 2.5, 3.5);
            let c = a - b;
            assert_eq!(c.px(), 0.5);
            assert_eq!(c.py(), 0.5);
            assert_eq!(c.pz(), 0.5);
            assert_eq!(c.e(), 0.5);
            assert_eq!(*c.kt2(), 0.5);
            assert_eq!(*c.phi(), 0.7853981633974483);
            assert_eq!(*c.rap(), 0.34657359027997264);
        }

        #[test]
        fn sub_assign_updates_components_and_reinitializes_cached_values() {
            let mut a = PseudoJet::new(1.0, 2.0, 3.0, 4.0);
            let b = PseudoJet::new(0.5, 1.5, 2.5, 3.5);
            a -= b;
            assert_eq!(a.px(), 0.5);
            assert_eq!(a.py(), 0.5);
            assert_eq!(a.pz(), 0.5);
            assert_eq!(a.e(), 0.5);
            assert_eq!(*a.kt2(), 0.5);
            assert_eq!(*a.phi(), 0.7853981633974483);
            assert_eq!(*a.rap(), 0.34657359027997264);
        }

        #[test]
        fn mul_scales_four_momentum_and_cached_kt2() {
            let jet: PseudoJet = PseudoJet::new(1.0, 2.0, 3.0, 4.0);
            let out: PseudoJet = jet * 3.0;
            assert_eq!(out.px(), 3.0);
            assert_eq!(out.py(), 6.0);
            assert_eq!(out.pz(), 9.0);
            assert_eq!(out.e(), 12.0);
            assert_eq!(*out.kt2(), 45.0);
            assert_eq!(*out.phi(), 1.1071487177940904);
            assert_eq!(*out.rap(), 0.9729550745276567);
        }

        #[test]
        fn mul_by_scalar_is_supported_from_left_side() {
            let jet = PseudoJet::new(1.0, 2.0, 3.0, 4.0);
            let out = 2.0 * jet;
            assert_eq!(out.px(), 2.0);
            assert_eq!(out.py(), 4.0);
            assert_eq!(out.pz(), 6.0);
            assert_eq!(out.e(), 8.0);
            assert_eq!(*out.kt2(), 20.0);
        }

        #[test]
        fn mul_assign_scales_in_place() {
            let mut jet = PseudoJet::new(1.0, 2.0, 3.0, 4.0);
            jet *= 2.0;
            assert_eq!(jet.px(), 2.0);
            assert_eq!(jet.py(), 4.0);
            assert_eq!(jet.pz(), 6.0);
            assert_eq!(jet.e(), 8.0);
            assert_eq!(*jet.kt2(), 20.0);
        }

        #[test]
        fn div_scales_by_inverse() {
            let jet = PseudoJet::new(2.0, 4.0, 6.0, 8.0);
            let out = jet / 2.0;
            assert_eq!(out.px(), 1.0);
            assert_eq!(out.py(), 2.0);
            assert_eq!(out.pz(), 3.0);
            assert_eq!(out.e(), 4.0);
            assert_eq!(*out.kt2(), 5.0);
        }

        #[test]
        fn div_assign_scales_by_inverse_in_place() {
            let mut jet = PseudoJet::new(2.0, 4.0, 6.0, 8.0);
            jet /= 2.0;
            assert_eq!(jet.px(), 1.0);
            assert_eq!(jet.py(), 2.0);
            assert_eq!(jet.pz(), 3.0);
            assert_eq!(jet.e(), 4.0);
            assert_eq!(*jet.kt2(), 5.0);
        }

        #[test]
        fn partial_eq_for_pseudojet_compares_four_momentum() {
            let a = PseudoJet::new(1.0, 2.0, 3.0, 4.0);
            let b = PseudoJet::new(1.0, 2.0, 3.0, 4.0);
            let c = PseudoJet::new(1.0, 2.0, 3.0, 5.0);
            assert!(a == b);
            assert!(a != c);
        }

        #[test]
        fn partial_eq_for_f64_supports_zero_only() {
            let zero = PseudoJet::new(0.0, 0.0, 0.0, 0.0);
            let non_zero = PseudoJet::new(1.0, 0.0, 0.0, 0.0);
            assert!(zero == 0.0);
            assert!(non_zero != 0.0);
        }

        #[test]
        #[should_panic(
            expected = "Comparing a PseudoJet with a non-zero constant (double) is not allowed"
        )]
        fn partial_eq_for_f64_panics_for_non_zero_constant() {
            let jet = PseudoJet::new(1.0, 2.0, 3.0, 4.0);
            let _ = jet == 1.0;
        }
    }
}
