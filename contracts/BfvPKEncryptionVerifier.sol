// SPDX-License-Identifier: MIT

pragma solidity ^0.8.0;

contract Halo2Verifier {
    uint256 internal constant    PROOF_LEN_CPTR = 0x44;
    uint256 internal constant        PROOF_CPTR = 0x64;
    uint256 internal constant NUM_INSTANCE_CPTR = 0x3284;
    uint256 internal constant     INSTANCE_CPTR = 0x32a4;

    uint256 internal constant FIRST_QUOTIENT_X_CPTR = 0x1124;
    uint256 internal constant  LAST_QUOTIENT_X_CPTR = 0x11a4;

    uint256 internal constant                VK_MPTR = 0x0d60;
    uint256 internal constant         VK_DIGEST_MPTR = 0x0d60;
    uint256 internal constant     NUM_INSTANCES_MPTR = 0x0d80;
    uint256 internal constant                 K_MPTR = 0x0da0;
    uint256 internal constant             N_INV_MPTR = 0x0dc0;
    uint256 internal constant             OMEGA_MPTR = 0x0de0;
    uint256 internal constant         OMEGA_INV_MPTR = 0x0e00;
    uint256 internal constant    OMEGA_INV_TO_L_MPTR = 0x0e20;
    uint256 internal constant   HAS_ACCUMULATOR_MPTR = 0x0e40;
    uint256 internal constant        ACC_OFFSET_MPTR = 0x0e60;
    uint256 internal constant     NUM_ACC_LIMBS_MPTR = 0x0e80;
    uint256 internal constant NUM_ACC_LIMB_BITS_MPTR = 0x0ea0;
    uint256 internal constant              G1_X_MPTR = 0x0ec0;
    uint256 internal constant              G1_Y_MPTR = 0x0ee0;
    uint256 internal constant            G2_X_1_MPTR = 0x0f00;
    uint256 internal constant            G2_X_2_MPTR = 0x0f20;
    uint256 internal constant            G2_Y_1_MPTR = 0x0f40;
    uint256 internal constant            G2_Y_2_MPTR = 0x0f60;
    uint256 internal constant      NEG_S_G2_X_1_MPTR = 0x0f80;
    uint256 internal constant      NEG_S_G2_X_2_MPTR = 0x0fa0;
    uint256 internal constant      NEG_S_G2_Y_1_MPTR = 0x0fc0;
    uint256 internal constant      NEG_S_G2_Y_2_MPTR = 0x0fe0;

    uint256 internal constant CHALLENGE_MPTR = 0x2040;

    uint256 internal constant THETA_MPTR = 0x2060;
    uint256 internal constant  BETA_MPTR = 0x2080;
    uint256 internal constant GAMMA_MPTR = 0x20a0;
    uint256 internal constant     Y_MPTR = 0x20c0;
    uint256 internal constant     X_MPTR = 0x20e0;
    uint256 internal constant  ZETA_MPTR = 0x2100;
    uint256 internal constant    NU_MPTR = 0x2120;
    uint256 internal constant    MU_MPTR = 0x2140;

    uint256 internal constant       ACC_LHS_X_MPTR = 0x2160;
    uint256 internal constant       ACC_LHS_Y_MPTR = 0x2180;
    uint256 internal constant       ACC_RHS_X_MPTR = 0x21a0;
    uint256 internal constant       ACC_RHS_Y_MPTR = 0x21c0;
    uint256 internal constant             X_N_MPTR = 0x21e0;
    uint256 internal constant X_N_MINUS_1_INV_MPTR = 0x2200;
    uint256 internal constant          L_LAST_MPTR = 0x2220;
    uint256 internal constant         L_BLIND_MPTR = 0x2240;
    uint256 internal constant             L_0_MPTR = 0x2260;
    uint256 internal constant   INSTANCE_EVAL_MPTR = 0x2280;
    uint256 internal constant   QUOTIENT_EVAL_MPTR = 0x22a0;
    uint256 internal constant      QUOTIENT_X_MPTR = 0x22c0;
    uint256 internal constant      QUOTIENT_Y_MPTR = 0x22e0;
    uint256 internal constant       G1_SCALAR_MPTR = 0x2300;
    uint256 internal constant   PAIRING_LHS_X_MPTR = 0x2320;
    uint256 internal constant   PAIRING_LHS_Y_MPTR = 0x2340;
    uint256 internal constant   PAIRING_RHS_X_MPTR = 0x2360;
    uint256 internal constant   PAIRING_RHS_Y_MPTR = 0x2380;

    function verifyProof(
        bytes calldata proof,
        uint256[] calldata instances
    ) public view returns (bool) {
        assembly {
            // Read EC point (x, y) at (proof_cptr, proof_cptr + 0x20),
            // and check if the point is on affine plane,
            // and store them in (hash_mptr, hash_mptr + 0x20).
            // Return updated (success, proof_cptr, hash_mptr).
            function read_ec_point(success, proof_cptr, hash_mptr, q) -> ret0, ret1, ret2 {
                let x := calldataload(proof_cptr)
                let y := calldataload(add(proof_cptr, 0x20))
                ret0 := and(success, lt(x, q))
                ret0 := and(ret0, lt(y, q))
                ret0 := and(ret0, eq(mulmod(y, y, q), addmod(mulmod(x, mulmod(x, x, q), q), 3, q)))
                mstore(hash_mptr, x)
                mstore(add(hash_mptr, 0x20), y)
                ret1 := add(proof_cptr, 0x40)
                ret2 := add(hash_mptr, 0x40)
            }

            // Squeeze challenge by keccak256(memory[0..hash_mptr]),
            // and store hash mod r as challenge in challenge_mptr,
            // and push back hash in 0x00 as the first input for next squeeze.
            // Return updated (challenge_mptr, hash_mptr).
            function squeeze_challenge(challenge_mptr, hash_mptr, r) -> ret0, ret1 {
                let hash := keccak256(0x00, hash_mptr)
                mstore(challenge_mptr, mod(hash, r))
                mstore(0x00, hash)
                ret0 := add(challenge_mptr, 0x20)
                ret1 := 0x20
            }

            // Squeeze challenge without absorbing new input from calldata,
            // by putting an extra 0x01 in memory[0x20] and squeeze by keccak256(memory[0..21]),
            // and store hash mod r as challenge in challenge_mptr,
            // and push back hash in 0x00 as the first input for next squeeze.
            // Return updated (challenge_mptr).
            function squeeze_challenge_cont(challenge_mptr, r) -> ret {
                mstore8(0x20, 0x01)
                let hash := keccak256(0x00, 0x21)
                mstore(challenge_mptr, mod(hash, r))
                mstore(0x00, hash)
                ret := add(challenge_mptr, 0x20)
            }

            // Batch invert values in memory[mptr_start..mptr_end] in place.
            // Return updated (success).
            function batch_invert(success, mptr_start, mptr_end, r) -> ret {
                let gp_mptr := mptr_end
                let gp := mload(mptr_start)
                let mptr := add(mptr_start, 0x20)
                for
                    {}
                    lt(mptr, sub(mptr_end, 0x20))
                    {}
                {
                    gp := mulmod(gp, mload(mptr), r)
                    mstore(gp_mptr, gp)
                    mptr := add(mptr, 0x20)
                    gp_mptr := add(gp_mptr, 0x20)
                }
                gp := mulmod(gp, mload(mptr), r)

                mstore(gp_mptr, 0x20)
                mstore(add(gp_mptr, 0x20), 0x20)
                mstore(add(gp_mptr, 0x40), 0x20)
                mstore(add(gp_mptr, 0x60), gp)
                mstore(add(gp_mptr, 0x80), sub(r, 2))
                mstore(add(gp_mptr, 0xa0), r)
                ret := and(success, staticcall(gas(), 0x05, gp_mptr, 0xc0, gp_mptr, 0x20))
                let all_inv := mload(gp_mptr)

                let first_mptr := mptr_start
                let second_mptr := add(first_mptr, 0x20)
                gp_mptr := sub(gp_mptr, 0x20)
                for
                    {}
                    lt(second_mptr, mptr)
                    {}
                {
                    let inv := mulmod(all_inv, mload(gp_mptr), r)
                    all_inv := mulmod(all_inv, mload(mptr), r)
                    mstore(mptr, inv)
                    mptr := sub(mptr, 0x20)
                    gp_mptr := sub(gp_mptr, 0x20)
                }
                let inv_first := mulmod(all_inv, mload(second_mptr), r)
                let inv_second := mulmod(all_inv, mload(first_mptr), r)
                mstore(first_mptr, inv_first)
                mstore(second_mptr, inv_second)
            }

            // Add (x, y) into point at (0x00, 0x20).
            // Return updated (success).
            function ec_add_acc(success, x, y) -> ret {
                mstore(0x40, x)
                mstore(0x60, y)
                ret := and(success, staticcall(gas(), 0x06, 0x00, 0x80, 0x00, 0x40))
            }

            // Scale point at (0x00, 0x20) by scalar.
            function ec_mul_acc(success, scalar) -> ret {
                mstore(0x40, scalar)
                ret := and(success, staticcall(gas(), 0x07, 0x00, 0x60, 0x00, 0x40))
            }

            // Add (x, y) into point at (0x80, 0xa0).
            // Return updated (success).
            function ec_add_tmp(success, x, y) -> ret {
                mstore(0xc0, x)
                mstore(0xe0, y)
                ret := and(success, staticcall(gas(), 0x06, 0x80, 0x80, 0x80, 0x40))
            }

            // Scale point at (0x80, 0xa0) by scalar.
            // Return updated (success).
            function ec_mul_tmp(success, scalar) -> ret {
                mstore(0xc0, scalar)
                ret := and(success, staticcall(gas(), 0x07, 0x80, 0x60, 0x80, 0x40))
            }

            // Perform pairing check.
            // Return updated (success).
            function ec_pairing(success, lhs_x, lhs_y, rhs_x, rhs_y) -> ret {
                mstore(0x00, lhs_x)
                mstore(0x20, lhs_y)
                mstore(0x40, mload(G2_X_1_MPTR))
                mstore(0x60, mload(G2_X_2_MPTR))
                mstore(0x80, mload(G2_Y_1_MPTR))
                mstore(0xa0, mload(G2_Y_2_MPTR))
                mstore(0xc0, rhs_x)
                mstore(0xe0, rhs_y)
                mstore(0x100, mload(NEG_S_G2_X_1_MPTR))
                mstore(0x120, mload(NEG_S_G2_X_2_MPTR))
                mstore(0x140, mload(NEG_S_G2_Y_1_MPTR))
                mstore(0x160, mload(NEG_S_G2_Y_2_MPTR))
                ret := and(success, staticcall(gas(), 0x08, 0x00, 0x180, 0x00, 0x20))
                ret := and(ret, mload(0x00))
            }

            // Modulus
            let q := 21888242871839275222246405745257275088696311157297823662689037894645226208583 // BN254 base field
            let r := 21888242871839275222246405745257275088548364400416034343698204186575808495617 // BN254 scalar field

            // Initialize success as true
            let success := true

            {
                // Load vk_digest and num_instances of vk into memory
                mstore(0x0d60, 0x220db9849286d79081212889b0d8f7c91016a6967be31ca804ab319d1916d8ca) // vk_digest
                mstore(0x0d80, 0x0000000000000000000000000000000000000000000000000000000000000004) // num_instances

                // Check valid length of proof
                success := and(success, eq(0x3220, calldataload(PROOF_LEN_CPTR)))

                // Check valid length of instances
                let num_instances := mload(NUM_INSTANCES_MPTR)
                success := and(success, eq(num_instances, calldataload(NUM_INSTANCE_CPTR)))

                // Absorb vk diegst
                mstore(0x00, mload(VK_DIGEST_MPTR))

                // Read instances and witness commitments and generate challenges
                let hash_mptr := 0x20
                let instance_cptr := INSTANCE_CPTR
                for
                    { let instance_cptr_end := add(instance_cptr, mul(0x20, num_instances)) }
                    lt(instance_cptr, instance_cptr_end)
                    {}
                {
                    let instance := calldataload(instance_cptr)
                    success := and(success, lt(instance, r))
                    mstore(hash_mptr, instance)
                    instance_cptr := add(instance_cptr, 0x20)
                    hash_mptr := add(hash_mptr, 0x20)
                }

                let proof_cptr := PROOF_CPTR
                let challenge_mptr := CHALLENGE_MPTR

                // Phase 1
                for
                    { let proof_cptr_end := add(proof_cptr, 0x80) }
                    lt(proof_cptr, proof_cptr_end)
                    {}
                {
                    success, proof_cptr, hash_mptr := read_ec_point(success, proof_cptr, hash_mptr, q)
                }

                challenge_mptr, hash_mptr := squeeze_challenge(challenge_mptr, hash_mptr, r)

                // Phase 2
                for
                    { let proof_cptr_end := add(proof_cptr, 0x07c0) }
                    lt(proof_cptr, proof_cptr_end)
                    {}
                {
                    success, proof_cptr, hash_mptr := read_ec_point(success, proof_cptr, hash_mptr, q)
                }

                challenge_mptr, hash_mptr := squeeze_challenge(challenge_mptr, hash_mptr, r)

                // Phase 3
                for
                    { let proof_cptr_end := add(proof_cptr, 0x0280) }
                    lt(proof_cptr, proof_cptr_end)
                    {}
                {
                    success, proof_cptr, hash_mptr := read_ec_point(success, proof_cptr, hash_mptr, q)
                }

                challenge_mptr, hash_mptr := squeeze_challenge(challenge_mptr, hash_mptr, r)
                challenge_mptr := squeeze_challenge_cont(challenge_mptr, r)

                // Phase 4
                for
                    { let proof_cptr_end := add(proof_cptr, 0x0600) }
                    lt(proof_cptr, proof_cptr_end)
                    {}
                {
                    success, proof_cptr, hash_mptr := read_ec_point(success, proof_cptr, hash_mptr, q)
                }

                challenge_mptr, hash_mptr := squeeze_challenge(challenge_mptr, hash_mptr, r)

                // Phase 5
                for
                    { let proof_cptr_end := add(proof_cptr, 0xc0) }
                    lt(proof_cptr, proof_cptr_end)
                    {}
                {
                    success, proof_cptr, hash_mptr := read_ec_point(success, proof_cptr, hash_mptr, q)
                }

                challenge_mptr, hash_mptr := squeeze_challenge(challenge_mptr, hash_mptr, r)

                // Read evaluations
                for
                    { let proof_cptr_end := add(proof_cptr, 0x2020) }
                    lt(proof_cptr, proof_cptr_end)
                    {}
                {
                    let eval := calldataload(proof_cptr)
                    success := and(success, lt(eval, r))
                    mstore(hash_mptr, eval)
                    proof_cptr := add(proof_cptr, 0x20)
                    hash_mptr := add(hash_mptr, 0x20)
                }

                // Read batch opening proof and generate challenges
                challenge_mptr, hash_mptr := squeeze_challenge(challenge_mptr, hash_mptr, r)       // zeta
                challenge_mptr := squeeze_challenge_cont(challenge_mptr, r)                        // nu

                success, proof_cptr, hash_mptr := read_ec_point(success, proof_cptr, hash_mptr, q) // W

                challenge_mptr, hash_mptr := squeeze_challenge(challenge_mptr, hash_mptr, r)       // mu

                success, proof_cptr, hash_mptr := read_ec_point(success, proof_cptr, hash_mptr, q) // W'

                // Load full vk into memory
                mstore(0x0d60, 0x220db9849286d79081212889b0d8f7c91016a6967be31ca804ab319d1916d8ca) // vk_digest
                mstore(0x0d80, 0x0000000000000000000000000000000000000000000000000000000000000004) // num_instances
                mstore(0x0da0, 0x000000000000000000000000000000000000000000000000000000000000000e) // k
                mstore(0x0dc0, 0x30638ce1a7661b6337a964756aa75257c6bf4778d89789ab819ce60c19b04001) // n_inv
                mstore(0x0de0, 0x2337acd19f40bf2b2aa212849e9a0c07d626d9ca335d73a09119dbe6eaab3cac) // omega
                mstore(0x0e00, 0x2f9c1d051b2a29bd1d13a09c1489aec5303c2fb2ac7d853ee7a58fdb65b90d7d) // omega_inv
                mstore(0x0e20, 0x190f2dc9aa5de00461a2889c7b743de1e467526db0fd664da4613c7ec12a9479) // omega_inv_to_l
                mstore(0x0e40, 0x0000000000000000000000000000000000000000000000000000000000000000) // has_accumulator
                mstore(0x0e60, 0x0000000000000000000000000000000000000000000000000000000000000000) // acc_offset
                mstore(0x0e80, 0x0000000000000000000000000000000000000000000000000000000000000000) // num_acc_limbs
                mstore(0x0ea0, 0x0000000000000000000000000000000000000000000000000000000000000000) // num_acc_limb_bits
                mstore(0x0ec0, 0x0000000000000000000000000000000000000000000000000000000000000001) // g1_x
                mstore(0x0ee0, 0x0000000000000000000000000000000000000000000000000000000000000002) // g1_y
                mstore(0x0f00, 0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2) // g2_x_1
                mstore(0x0f20, 0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed) // g2_x_2
                mstore(0x0f40, 0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b) // g2_y_1
                mstore(0x0f60, 0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa) // g2_y_2
                mstore(0x0f80, 0x0181624e80f3d6ae28df7e01eaeab1c0e919877a3b8a6b7fbc69a6817d596ea2) // neg_s_g2_x_1
                mstore(0x0fa0, 0x1783d30dcb12d259bb89098addf6280fa4b653be7a152542a28f7b926e27e648) // neg_s_g2_x_2
                mstore(0x0fc0, 0x00ae44489d41a0d179e2dfdc03bddd883b7109f8b6ae316a59e815c1a6b35304) // neg_s_g2_y_1
                mstore(0x0fe0, 0x0b2147ab62a386bd63e6de1522109b8c9588ab466f5aadfde8c41ca3749423ee) // neg_s_g2_y_2
                mstore(0x1000, 0x2b5ea13bc24e2b5ecc529ee04eb630aa7c9336a149e7124c1d1642da8a70bec3) // fixed_comms[0].x
                mstore(0x1020, 0x03d04df3861d2ea92c7318a09463f9ff4b169ffe2933ee0deee139a519530983) // fixed_comms[0].y
                mstore(0x1040, 0x04d0e86e8aa5612353d11f0f5877f0a8ae07dbe3cb0fb046b3262820f1f79495) // fixed_comms[1].x
                mstore(0x1060, 0x14516cf1b41f08582e3127c7ff8fcef84ca35af040cc1b5b2bf3e30e4a81e4c8) // fixed_comms[1].y
                mstore(0x1080, 0x0000000000000000000000000000000000000000000000000000000000000000) // fixed_comms[2].x
                mstore(0x10a0, 0x0000000000000000000000000000000000000000000000000000000000000000) // fixed_comms[2].y
                mstore(0x10c0, 0x0000000000000000000000000000000000000000000000000000000000000000) // fixed_comms[3].x
                mstore(0x10e0, 0x0000000000000000000000000000000000000000000000000000000000000000) // fixed_comms[3].y
                mstore(0x1100, 0x17fa79a50f3c69260545833791570f97d564a57894ac289f684f582e9d1681ec) // fixed_comms[4].x
                mstore(0x1120, 0x1be4fb19bf76593d235d0e469c680bdffa22afc54a786acb7179caaa32e2eb76) // fixed_comms[4].y
                mstore(0x1140, 0x1fa62e48b2d88aca2a35e47e9d0e46eb71551015edd04591e30a29fa88635177) // fixed_comms[5].x
                mstore(0x1160, 0x0b4ca55515949e844ea8fe91f1b5ef270f2528bbb310f4d14beb4ba16a2a7300) // fixed_comms[5].y
                mstore(0x1180, 0x0cc8cc5711956c11761195d4cbec13fa292294cddcfc12bf93b21770d5d612b4) // fixed_comms[6].x
                mstore(0x11a0, 0x2544135a54fbe01df6b59c01734370701691d2a60eaa9af485e837c217bbaf00) // fixed_comms[6].y
                mstore(0x11c0, 0x03697ada88bec0192e78d06025413014830d628f8e1bbf99c3f19c536915e54e) // fixed_comms[7].x
                mstore(0x11e0, 0x01d3bce2cca3657327e13cd8f97f4ad94e459bd60d24e01bd68280d291cd0889) // fixed_comms[7].y
                mstore(0x1200, 0x23a905a9f1a02421776b2958a56bbfcfdf82d840e8da3be15de647e72ff21005) // fixed_comms[8].x
                mstore(0x1220, 0x175eca1100d1a0edb4e0dc26b12cd4143b6ff9db47cf36e1882e5752d162acf5) // fixed_comms[8].y
                mstore(0x1240, 0x25c7f9d60a9d01eba375a9296c1331e6f89a7b6616c01f2d067f0545e96a140e) // fixed_comms[9].x
                mstore(0x1260, 0x18e7db4f211befccad649e3ec6c1ce064e8c16cefeb69390a93e9ca8cc5e40f4) // fixed_comms[9].y
                mstore(0x1280, 0x03fa467b5cbf3fb0c2092a08528aeb8fbaf3db64172d86397603e84e5e2cf253) // fixed_comms[10].x
                mstore(0x12a0, 0x259a7832d4f4bf8fd16e978b74c688adee116502bbacdd294dc4d28b6450e3c6) // fixed_comms[10].y
                mstore(0x12c0, 0x0f30107f996c8dd3d166b0364963b1c9d1e8e5ad82f6fc080243120211130c8a) // fixed_comms[11].x
                mstore(0x12e0, 0x0924166e0da0933439a1fa970ee50483ae9d69fc25ce6d964d04208dd5b58cfa) // fixed_comms[11].y
                mstore(0x1300, 0x15c3671f4750e6089b8c3c4bb778f23aa1976f14ed8d222dd8b6d6bef7ae838a) // fixed_comms[12].x
                mstore(0x1320, 0x231ab059df52977aa2c56f3a217f899033a11bf333fa45f745a238215968bd47) // fixed_comms[12].y
                mstore(0x1340, 0x2996c1a53e67945543bd51f28aec091547aeb87891811aa310d9215705214dcf) // fixed_comms[13].x
                mstore(0x1360, 0x2e853300a34d0394609b59132bd437ac2bb094715e558b75a21f793df3f793a8) // fixed_comms[13].y
                mstore(0x1380, 0x3041b8018ed2e841d5725ed18690751e14329e8a86a249c22583e5a4140fb69b) // fixed_comms[14].x
                mstore(0x13a0, 0x09fa4ac9d95132cbb7a3fed063b736bc3f1801de5048d1a0ae3d4b72d50d2748) // fixed_comms[14].y
                mstore(0x13c0, 0x2e6e539f8f516527f2caccd7f782fbf4b88e7c565b036262d3a5456663835143) // fixed_comms[15].x
                mstore(0x13e0, 0x07c1c3860b76d683b19af0066687683ad3e5fd6d107e67f7d3d65317d0eab543) // fixed_comms[15].y
                mstore(0x1400, 0x2d293cf80417ca9bcb9404c5b3aa8f71c684c5946dc4ca29ac53ec6bc82fdb7d) // fixed_comms[16].x
                mstore(0x1420, 0x29cb083173c1b0cde972c20cd202573ff1f5266642ac8335561410486be72b62) // fixed_comms[16].y
                mstore(0x1440, 0x2f34ee4837045d332a8054d6c43d0180c4e60b14c8e02f6c2bb12ef6777df748) // fixed_comms[17].x
                mstore(0x1460, 0x12a73aa31f456e5ee0df93dc1ca2a82def82a2661da5fa02a1b8317078a16ebe) // fixed_comms[17].y
                mstore(0x1480, 0x0ba670950511a895b66086ff5796a06ec3561d854dd3b1645f5dd5bf61f1bda6) // fixed_comms[18].x
                mstore(0x14a0, 0x1532637ccdafb1da847c595c7612e591e02d2e80602ab081e63579138087e13e) // fixed_comms[18].y
                mstore(0x14c0, 0x1e78e6299d1f4088b50df22666327fb5f251dd79930f850ba958993de20c5dd0) // fixed_comms[19].x
                mstore(0x14e0, 0x13abcc8bf711259399409d6a6490de7eacd1f938f7da60dd840711723d7ac889) // fixed_comms[19].y
                mstore(0x1500, 0x18810c315ba9c59d9a10b1cd3035f56ed57a08b98ed026fd45230de149adf9b8) // fixed_comms[20].x
                mstore(0x1520, 0x21582fddd8c9fa81f336138010c9aa5e86a062e2c4c6e7cca44cba10ae98aaa0) // fixed_comms[20].y
                mstore(0x1540, 0x06b365bd2ff8febd1c44fadcaf99d233f9ad1d7b2c846f92855b7151653d960c) // fixed_comms[21].x
                mstore(0x1560, 0x02a2bd0251e33d7114ccf0c4f60c5e3311c70c8b679b392436c8edadbf1c6615) // fixed_comms[21].y
                mstore(0x1580, 0x1ac146b095610f90e12674df5399e007dfe010c37015725f373c5d2b03ea5744) // fixed_comms[22].x
                mstore(0x15a0, 0x10f1cca6ff9224fe51d0a28f6907a8bffde2bf0665c4cb123fbfd9006b51807e) // fixed_comms[22].y
                mstore(0x15c0, 0x1d7468bd368ba342a47cd0ae84014d10b2a8afea8efff2ef74c7d96502f7d744) // fixed_comms[23].x
                mstore(0x15e0, 0x06d7380bf09b863b7c5f194943142695e8051aa5126e76fd593793776576a821) // fixed_comms[23].y
                mstore(0x1600, 0x0cc8cc5711956c11761195d4cbec13fa292294cddcfc12bf93b21770d5d612b4) // fixed_comms[24].x
                mstore(0x1620, 0x2544135a54fbe01df6b59c01734370701691d2a60eaa9af485e837c217bbaf00) // fixed_comms[24].y
                mstore(0x1640, 0x2e0052f46b21287fb187fecfcb6a4da84f41192bb3691082efcf3aa0aca5d7ee) // fixed_comms[25].x
                mstore(0x1660, 0x1fdf1913c2d610bacd971eddad9ea23f8c0c5633dfb5ea5f292cc8670dce5f2e) // fixed_comms[25].y
                mstore(0x1680, 0x18a62fc66205835dde6bd76883ff2fc576734070d0d7610c779a877ab4b3d634) // fixed_comms[26].x
                mstore(0x16a0, 0x0772ac519c5a8dbc135394d1c1e9b4fd863e57a8691a0f5fe7d5af041fbde44d) // fixed_comms[26].y
                mstore(0x16c0, 0x16264ec60aa8d19ed483f70736ccad9d8f09fd2a399832bd229de10a728c7bfc) // fixed_comms[27].x
                mstore(0x16e0, 0x04cc4535a3251dfb716a49ef1b995a13997e8251e87251ed7c4e16c97cbce645) // fixed_comms[27].y
                mstore(0x1700, 0x116de0c5670c48600a9d06820999c26f7820b4831e61cbb2996e497a07042d7c) // fixed_comms[28].x
                mstore(0x1720, 0x1dc2014cb6fc7088b95c15f6dcb0d2d8e9169c6d304adac48e5507f2a2fc3371) // fixed_comms[28].y
                mstore(0x1740, 0x04606af192bd8526e3583b20e5d3408c19ae64f410c59f33d62d411fbe63f278) // fixed_comms[29].x
                mstore(0x1760, 0x13458f142c14515eb3f0636e2c1dbc556092f135caaa33cf78da74a5e203ff27) // fixed_comms[29].y
                mstore(0x1780, 0x1fdb5c5225c5afe2629e1bf89d95b01dd85358a96eb65bccb4598fda14477536) // permutation_comms[0].x
                mstore(0x17a0, 0x2fe0b6f788ba24a51ebd6a538ffd37dee471057a1cf1e34fc14c4cb912f09b65) // permutation_comms[0].y
                mstore(0x17c0, 0x14f942ffdce53ae8484ef0d70a6c577130c9bafd1e3780c7be76b14a34acae84) // permutation_comms[1].x
                mstore(0x17e0, 0x0ea82b2658fd277d5ae20c9473aa5827e55abf69af0b7acabb7ceae7105bd5fc) // permutation_comms[1].y
                mstore(0x1800, 0x09ef66b16f8c2a205aa4f9ed75ecb63394d6a04bfbcd58b0b0cfda7d6b9f9f33) // permutation_comms[2].x
                mstore(0x1820, 0x249f6c742bb596f8d5b69e8f48c384b165a1ace3fe38fc27428f117a44bc33d8) // permutation_comms[2].y
                mstore(0x1840, 0x20fda67fa2e86797a45ab06cd26e2bf4ae12ee0685add7c215f2eab15e2842b4) // permutation_comms[3].x
                mstore(0x1860, 0x256a32b6ce557e285d5d9bd3454d923a869f3eaf38ffa55bf8d335d6a77f3408) // permutation_comms[3].y
                mstore(0x1880, 0x0cdc044494bc84973c3a9dc749c5e4cec5b261ecfa5af666b86712b5d20469ed) // permutation_comms[4].x
                mstore(0x18a0, 0x16c3b22690046dcc1311d09e224a299bc2d6c40c0723b648a4687641fdcef2cf) // permutation_comms[4].y
                mstore(0x18c0, 0x239a61f1ee1f39cd1415cef692f43551ec300dd4eed70f044e116082262f133b) // permutation_comms[5].x
                mstore(0x18e0, 0x25842a7f1aad9f1b42761db3a53c084cbfcd6d5ad457d9656b1b37c9d992a115) // permutation_comms[5].y
                mstore(0x1900, 0x2f341cafe2f8f6f5933e845cf2229435531b4fd0527879a3ad32784cee8846a3) // permutation_comms[6].x
                mstore(0x1920, 0x1f379dbecc00f3a99914c87b2f88eb6c34615ef6d8bf5e06a64deafea543ded1) // permutation_comms[6].y
                mstore(0x1940, 0x0f918c4f5e48130d8073bbb2ab7c99e788e96b401cfda62e3b9312490adc411f) // permutation_comms[7].x
                mstore(0x1960, 0x1d8a91c8f60522a4011c8401777b9df88e8b191ae8138a5a1c4d7b2059c215bb) // permutation_comms[7].y
                mstore(0x1980, 0x20e65218df9fb7faa68c23ab0476b233635d82c4d3c46e2ebbe7cb004d91552d) // permutation_comms[8].x
                mstore(0x19a0, 0x2b076c56f26090319f757972352ceb03e41fa718e3c619357d35f5448332889a) // permutation_comms[8].y
                mstore(0x19c0, 0x168e53e3be500ef8b4c9cf90a5c1cf422f101a973a338e64d860d7cfbe5a79f3) // permutation_comms[9].x
                mstore(0x19e0, 0x0a844d4b912772e8e6fdf3b490702441bf84ace042a145e208114dfb801e4e1b) // permutation_comms[9].y
                mstore(0x1a00, 0x1578d40ef72e444189101fe13c615b63ada531253ea7eaa55f61dd742568175a) // permutation_comms[10].x
                mstore(0x1a20, 0x053860610064c2d77800fc9ba89887783d53af2df36a7ce18f0e9f749451a83e) // permutation_comms[10].y
                mstore(0x1a40, 0x1d8e1a34c95039eb14366e6ebf4e2df121b661946ad268f52a6ff573d8013bb4) // permutation_comms[11].x
                mstore(0x1a60, 0x107f32928d662ede1f4fb1ccfcb915e651fd85cd8112b703332dc19e8145fe4a) // permutation_comms[11].y
                mstore(0x1a80, 0x0196e34e83a750ea86b70562f252098f905438d81d705d3a8c1980dfb43ccc33) // permutation_comms[12].x
                mstore(0x1aa0, 0x0febc098e0a49a491a43096abae01a9851b7b398102181aa82a957c7b1f1685a) // permutation_comms[12].y
                mstore(0x1ac0, 0x131eb8c708f24e37d64a952ce8fb3ceb3efcabebf07453dcb676ff43b085c2f3) // permutation_comms[13].x
                mstore(0x1ae0, 0x1166664d30857d79de9a858e8ff7bbaa9294145bfcf38866475793d518f98bb4) // permutation_comms[13].y
                mstore(0x1b00, 0x27c5090cf56aab972eba3645e4d24b4ec5906ec6ed62f1a376b09d355089292d) // permutation_comms[14].x
                mstore(0x1b20, 0x0890d798b56cab6c86869d534031367fcf7144f64d2210be7771cbea86618c08) // permutation_comms[14].y
                mstore(0x1b40, 0x282da800056665cdcf6128c23d0e6010b3e70c89ec75324f095643aee600acb8) // permutation_comms[15].x
                mstore(0x1b60, 0x006229f92eeaa280360331c6f1346514e0fb1c15ec41ff15f87cc25de65df7a8) // permutation_comms[15].y
                mstore(0x1b80, 0x05b6e335c3a546d22ccf0c57bf17424fc0ff5ec90acb951b4821f2ab40bf23a3) // permutation_comms[16].x
                mstore(0x1ba0, 0x2dea244c41486536c85b06b614a96560375b5463b7ca0e422d0dd6205c4e16b7) // permutation_comms[16].y
                mstore(0x1bc0, 0x17d29e29a733bfaf11c90ed82529d3e41077278cc6b498acddfaacf4746a2fa4) // permutation_comms[17].x
                mstore(0x1be0, 0x1ce5a212210ac2ef413aa8ea6793422972c148ca411564e815c8d3fa0ceb42cb) // permutation_comms[17].y
                mstore(0x1c00, 0x1b7452587ca72de16fb073916703b84c30bee6f4d2b6da56f3583dab5c09899d) // permutation_comms[18].x
                mstore(0x1c20, 0x2818ec531c289b83708bd8ec01d82556376ef74abc635dbc6eda20b9aa98396c) // permutation_comms[18].y
                mstore(0x1c40, 0x1527482fca8c3dafd53dc5081fc0f98c42c81b5f3857ce977f2d672d503f0f69) // permutation_comms[19].x
                mstore(0x1c60, 0x02ecbe5477eb4c3fc292812cef5b9bdb126394a2c22b39f795fe2f70b525e9a7) // permutation_comms[19].y
                mstore(0x1c80, 0x004c6386e5863f9c1665d1a153bc7e98808c2c71f8b715b2c12962e98e97bd6a) // permutation_comms[20].x
                mstore(0x1ca0, 0x1aa93abdb410aa01316c4d4a5612f78c3d22eb19f63604c4157a7d9d30fdae34) // permutation_comms[20].y
                mstore(0x1cc0, 0x15e8325b18f208ccef39c9698c9f2d557eef5027531271223ee0fe6c79b78f5c) // permutation_comms[21].x
                mstore(0x1ce0, 0x0b3ac1486107fe4fda78533256a3572cd151529d31e2131ca897201b0db66b26) // permutation_comms[21].y
                mstore(0x1d00, 0x1e7292703c6cb3a9e316ab454f75ee8a00a3e630f765975860d87a8029a61494) // permutation_comms[22].x
                mstore(0x1d20, 0x060e41d2d9ec61b2f101c777151cc11cede26f7d44d739abd4d703e2132f278f) // permutation_comms[22].y
                mstore(0x1d40, 0x15163c37c7091ef99c9d1951a57da2151ede1becc25ce136acd6d676c010413b) // permutation_comms[23].x
                mstore(0x1d60, 0x191bfdeb08e3110922042a7075ce0e1b2c014124b1ca68cf62507b301e399408) // permutation_comms[23].y
                mstore(0x1d80, 0x295e3e088aa7649b613ca44922cee28fe3e3243661a4f334b250c3479b7dd7fa) // permutation_comms[24].x
                mstore(0x1da0, 0x26fc25b160757454ed6190ec463b825d23a8d2b32a0c98cc9329fba55dd1881b) // permutation_comms[24].y
                mstore(0x1dc0, 0x28872ccba33d8c439db08a8b27d9ae57678d4a2e239fbb18d53b0fb91318e075) // permutation_comms[25].x
                mstore(0x1de0, 0x0eda18d077b2ce0051ec2e41b3e277e5ba963d96913c024624b6a54fdb15306f) // permutation_comms[25].y
                mstore(0x1e00, 0x18767a85f8f9a15eb0af8976f1daaedcdd9a44647d17bec456113af25fd197e2) // permutation_comms[26].x
                mstore(0x1e20, 0x23a770c5b9250bf951fd27bf9cb27fe227c62f481c4415311b14bb4e3aaa346d) // permutation_comms[26].y
                mstore(0x1e40, 0x2f36acb39a0eafe3c2ede7e918ecd2a564e4f1b2a60f77b1ac0e802775a45e7f) // permutation_comms[27].x
                mstore(0x1e60, 0x2d723a9e5d714d949307c328daebab919e6b96436149353f79c4d1afd61214af) // permutation_comms[27].y
                mstore(0x1e80, 0x1a1f7ba1e2c83f477054eb20f885aa84f05ef0d6894f3f0a119a7354c7cfcae5) // permutation_comms[28].x
                mstore(0x1ea0, 0x01bf01604b8d181e11520345304353aa9ae087869f011999decdca6f44ad1906) // permutation_comms[28].y
                mstore(0x1ec0, 0x184f31437e2893da907d9a3c91db13f394d28b684cf3995c75cb4e82d0d35630) // permutation_comms[29].x
                mstore(0x1ee0, 0x1f65cd4724252e16314b8908bdd13f9425d3fed82782df4fce72efa9b2f315d6) // permutation_comms[29].y
                mstore(0x1f00, 0x2e466e6da144644c6ad46a564bebc79d1863600c2681777f73f01a1d36dfff04) // permutation_comms[30].x
                mstore(0x1f20, 0x050db1efe615ec22b065c263772844d5e137fb36f325a766b1eb00ec143415ff) // permutation_comms[30].y
                mstore(0x1f40, 0x1d891e530368ccdc261659d104c06a7ebdfe5cec4067af776bb8b5f633752415) // permutation_comms[31].x
                mstore(0x1f60, 0x304cc61627e3eb9475dfe5db0e2d470f4d641257dae0be2b38cbeaa84414155e) // permutation_comms[31].y
                mstore(0x1f80, 0x2733b5b632a8396880f16d9bfe3d4e035a5460a4f5e4b1e4fbb42b462d9d1238) // permutation_comms[32].x
                mstore(0x1fa0, 0x2b48a86a2d61b42dec94bfee0361c2584f2ba006daeac06460eeabb6c881ae4d) // permutation_comms[32].y
                mstore(0x1fc0, 0x25d7c9b7bdfbdcd856d0577389d48f40afca716314f6c8299070372cbd8e7fba) // permutation_comms[33].x
                mstore(0x1fe0, 0x1a57686bab77fd9bb1b81147a70c919ed186e6b2bb0e99f1506c4677486b7108) // permutation_comms[33].y
                mstore(0x2000, 0x04e103247a214f69be096012ea9bcfd0dfe0b2aef0e7b4748a1efb22e42f0ec9) // permutation_comms[34].x
                mstore(0x2020, 0x249e0895536da32b194cd06e02eed87f7ec4740b1a526828d6a135fbeca48baf) // permutation_comms[34].y

                // Read accumulator from instances
                if mload(HAS_ACCUMULATOR_MPTR) {
                    let num_limbs := mload(NUM_ACC_LIMBS_MPTR)
                    let num_limb_bits := mload(NUM_ACC_LIMB_BITS_MPTR)

                    let cptr := add(INSTANCE_CPTR, mul(mload(ACC_OFFSET_MPTR), 0x20))
                    let lhs_y_off := mul(num_limbs, 0x20)
                    let rhs_x_off := mul(lhs_y_off, 2)
                    let rhs_y_off := mul(lhs_y_off, 3)
                    let lhs_x := calldataload(cptr)
                    let lhs_y := calldataload(add(cptr, lhs_y_off))
                    let rhs_x := calldataload(add(cptr, rhs_x_off))
                    let rhs_y := calldataload(add(cptr, rhs_y_off))
                    for
                        {
                            let cptr_end := add(cptr, mul(0x20, num_limbs))
                            let shift := num_limb_bits
                        }
                        lt(cptr, cptr_end)
                        {}
                    {
                        cptr := add(cptr, 0x20)
                        lhs_x := add(lhs_x, shl(shift, calldataload(cptr)))
                        lhs_y := add(lhs_y, shl(shift, calldataload(add(cptr, lhs_y_off))))
                        rhs_x := add(rhs_x, shl(shift, calldataload(add(cptr, rhs_x_off))))
                        rhs_y := add(rhs_y, shl(shift, calldataload(add(cptr, rhs_y_off))))
                        shift := add(shift, num_limb_bits)
                    }

                    success := and(success, and(lt(lhs_x, q), lt(lhs_y, q)))
                    success := and(success, eq(mulmod(lhs_y, lhs_y, q), addmod(mulmod(lhs_x, mulmod(lhs_x, lhs_x, q), q), 3, q)))
                    success := and(success, and(lt(rhs_x, q), lt(rhs_y, q)))
                    success := and(success, eq(mulmod(rhs_y, rhs_y, q), addmod(mulmod(rhs_x, mulmod(rhs_x, rhs_x, q), q), 3, q)))

                    mstore(ACC_LHS_X_MPTR, lhs_x)
                    mstore(ACC_LHS_Y_MPTR, lhs_y)
                    mstore(ACC_RHS_X_MPTR, rhs_x)
                    mstore(ACC_RHS_Y_MPTR, rhs_y)
                }

                pop(q)
            }

            // Revert earlier if anything from calldata is invalid
            if iszero(success) {
                revert(0, 0)
            }

            // Compute lagrange evaluations and instance evaluation
            {
                let k := mload(K_MPTR)
                let x := mload(X_MPTR)
                let x_n := x
                for
                    { let idx := 0 }
                    lt(idx, k)
                    { idx := add(idx, 1) }
                {
                    x_n := mulmod(x_n, x_n, r)
                }

                let omega := mload(OMEGA_MPTR)

                let mptr := X_N_MPTR
                let mptr_end := add(mptr, mul(0x20, add(mload(NUM_INSTANCES_MPTR), 7)))
                if iszero(mload(NUM_INSTANCES_MPTR)) {
                    mptr_end := add(mptr_end, 0x20)
                }
                for
                    { let pow_of_omega := mload(OMEGA_INV_TO_L_MPTR) }
                    lt(mptr, mptr_end)
                    { mptr := add(mptr, 0x20) }
                {
                    mstore(mptr, addmod(x, sub(r, pow_of_omega), r))
                    pow_of_omega := mulmod(pow_of_omega, omega, r)
                }
                let x_n_minus_1 := addmod(x_n, sub(r, 1), r)
                mstore(mptr_end, x_n_minus_1)
                success := batch_invert(success, X_N_MPTR, add(mptr_end, 0x20), r)

                mptr := X_N_MPTR
                let l_i_common := mulmod(x_n_minus_1, mload(N_INV_MPTR), r)
                for
                    { let pow_of_omega := mload(OMEGA_INV_TO_L_MPTR) }
                    lt(mptr, mptr_end)
                    { mptr := add(mptr, 0x20) }
                {
                    mstore(mptr, mulmod(l_i_common, mulmod(mload(mptr), pow_of_omega, r), r))
                    pow_of_omega := mulmod(pow_of_omega, omega, r)
                }

                let l_blind := mload(add(X_N_MPTR, 0x20))
                let l_i_cptr := add(X_N_MPTR, 0x40)
                for
                    { let l_i_cptr_end := add(X_N_MPTR, 0xe0) }
                    lt(l_i_cptr, l_i_cptr_end)
                    { l_i_cptr := add(l_i_cptr, 0x20) }
                {
                    l_blind := addmod(l_blind, mload(l_i_cptr), r)
                }

                let instance_eval := 0
                for
                    {
                        let instance_cptr := INSTANCE_CPTR
                        let instance_cptr_end := add(instance_cptr, mul(0x20, mload(NUM_INSTANCES_MPTR)))
                    }
                    lt(instance_cptr, instance_cptr_end)
                    {
                        instance_cptr := add(instance_cptr, 0x20)
                        l_i_cptr := add(l_i_cptr, 0x20)
                    }
                {
                    instance_eval := addmod(instance_eval, mulmod(mload(l_i_cptr), calldataload(instance_cptr), r), r)
                }

                let x_n_minus_1_inv := mload(mptr_end)
                let l_last := mload(X_N_MPTR)
                let l_0 := mload(add(X_N_MPTR, 0xe0))

                mstore(X_N_MPTR, x_n)
                mstore(X_N_MINUS_1_INV_MPTR, x_n_minus_1_inv)
                mstore(L_LAST_MPTR, l_last)
                mstore(L_BLIND_MPTR, l_blind)
                mstore(L_0_MPTR, l_0)
                mstore(INSTANCE_EVAL_MPTR, instance_eval)
            }

            // Compute quotient evavluation
            {
                let quotient_eval_numer
                let delta := 4131629893567559867359510883348571134090853742863529169391034518566172092834
                let y := mload(Y_MPTR)
                {
                    let f_2 := calldataload(0x2044)
                    let a_0 := calldataload(0x11e4)
                    let a_0_next_1 := calldataload(0x1204)
                    let a_0_next_2 := calldataload(0x1224)
                    let var0 := mulmod(a_0_next_1, a_0_next_2, r)
                    let var1 := addmod(a_0, var0, r)
                    let a_0_next_3 := calldataload(0x1244)
                    let var2 := sub(r, a_0_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_2, var3, r)
                    quotient_eval_numer := var4
                }
                {
                    let f_3 := calldataload(0x2064)
                    let a_1 := calldataload(0x1264)
                    let a_1_next_1 := calldataload(0x1284)
                    let a_1_next_2 := calldataload(0x12a4)
                    let var0 := mulmod(a_1_next_1, a_1_next_2, r)
                    let var1 := addmod(a_1, var0, r)
                    let a_1_next_3 := calldataload(0x12c4)
                    let var2 := sub(r, a_1_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_3, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_4 := calldataload(0x2084)
                    let a_2 := calldataload(0x12e4)
                    let a_2_next_1 := calldataload(0x1304)
                    let a_2_next_2 := calldataload(0x1324)
                    let var0 := mulmod(a_2_next_1, a_2_next_2, r)
                    let var1 := addmod(a_2, var0, r)
                    let a_2_next_3 := calldataload(0x1344)
                    let var2 := sub(r, a_2_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_4, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_5 := calldataload(0x20a4)
                    let a_3 := calldataload(0x1364)
                    let a_3_next_1 := calldataload(0x1384)
                    let a_3_next_2 := calldataload(0x13a4)
                    let var0 := mulmod(a_3_next_1, a_3_next_2, r)
                    let var1 := addmod(a_3, var0, r)
                    let a_3_next_3 := calldataload(0x13c4)
                    let var2 := sub(r, a_3_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_5, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_6 := calldataload(0x20c4)
                    let a_4 := calldataload(0x13e4)
                    let a_4_next_1 := calldataload(0x1404)
                    let a_4_next_2 := calldataload(0x1424)
                    let var0 := mulmod(a_4_next_1, a_4_next_2, r)
                    let var1 := addmod(a_4, var0, r)
                    let a_4_next_3 := calldataload(0x1444)
                    let var2 := sub(r, a_4_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_6, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_7 := calldataload(0x20e4)
                    let a_5 := calldataload(0x1464)
                    let a_5_next_1 := calldataload(0x1484)
                    let a_5_next_2 := calldataload(0x14a4)
                    let var0 := mulmod(a_5_next_1, a_5_next_2, r)
                    let var1 := addmod(a_5, var0, r)
                    let a_5_next_3 := calldataload(0x14c4)
                    let var2 := sub(r, a_5_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_7, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_8 := calldataload(0x2104)
                    let a_6 := calldataload(0x14e4)
                    let a_6_next_1 := calldataload(0x1504)
                    let a_6_next_2 := calldataload(0x1524)
                    let var0 := mulmod(a_6_next_1, a_6_next_2, r)
                    let var1 := addmod(a_6, var0, r)
                    let a_6_next_3 := calldataload(0x1544)
                    let var2 := sub(r, a_6_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_8, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_9 := calldataload(0x2124)
                    let a_7 := calldataload(0x1564)
                    let a_7_next_1 := calldataload(0x1584)
                    let a_7_next_2 := calldataload(0x15a4)
                    let var0 := mulmod(a_7_next_1, a_7_next_2, r)
                    let var1 := addmod(a_7, var0, r)
                    let a_7_next_3 := calldataload(0x15c4)
                    let var2 := sub(r, a_7_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_9, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_10 := calldataload(0x2144)
                    let a_8 := calldataload(0x15e4)
                    let a_8_next_1 := calldataload(0x1604)
                    let a_8_next_2 := calldataload(0x1624)
                    let var0 := mulmod(a_8_next_1, a_8_next_2, r)
                    let var1 := addmod(a_8, var0, r)
                    let a_8_next_3 := calldataload(0x1644)
                    let var2 := sub(r, a_8_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_10, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_11 := calldataload(0x2164)
                    let a_9 := calldataload(0x1664)
                    let a_9_next_1 := calldataload(0x1684)
                    let a_9_next_2 := calldataload(0x16a4)
                    let var0 := mulmod(a_9_next_1, a_9_next_2, r)
                    let var1 := addmod(a_9, var0, r)
                    let a_9_next_3 := calldataload(0x16c4)
                    let var2 := sub(r, a_9_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_11, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_12 := calldataload(0x2184)
                    let a_10 := calldataload(0x16e4)
                    let a_10_next_1 := calldataload(0x1704)
                    let a_10_next_2 := calldataload(0x1724)
                    let var0 := mulmod(a_10_next_1, a_10_next_2, r)
                    let var1 := addmod(a_10, var0, r)
                    let a_10_next_3 := calldataload(0x1744)
                    let var2 := sub(r, a_10_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_12, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_13 := calldataload(0x21a4)
                    let a_11 := calldataload(0x1764)
                    let a_11_next_1 := calldataload(0x1784)
                    let a_11_next_2 := calldataload(0x17a4)
                    let var0 := mulmod(a_11_next_1, a_11_next_2, r)
                    let var1 := addmod(a_11, var0, r)
                    let a_11_next_3 := calldataload(0x17c4)
                    let var2 := sub(r, a_11_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_13, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_14 := calldataload(0x21c4)
                    let a_12 := calldataload(0x17e4)
                    let a_12_next_1 := calldataload(0x1804)
                    let a_12_next_2 := calldataload(0x1824)
                    let var0 := mulmod(a_12_next_1, a_12_next_2, r)
                    let var1 := addmod(a_12, var0, r)
                    let a_12_next_3 := calldataload(0x1844)
                    let var2 := sub(r, a_12_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_14, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_15 := calldataload(0x21e4)
                    let a_13 := calldataload(0x1864)
                    let a_13_next_1 := calldataload(0x1884)
                    let a_13_next_2 := calldataload(0x18a4)
                    let var0 := mulmod(a_13_next_1, a_13_next_2, r)
                    let var1 := addmod(a_13, var0, r)
                    let a_13_next_3 := calldataload(0x18c4)
                    let var2 := sub(r, a_13_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_15, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_16 := calldataload(0x2204)
                    let a_14 := calldataload(0x18e4)
                    let a_14_next_1 := calldataload(0x1904)
                    let a_14_next_2 := calldataload(0x1924)
                    let var0 := mulmod(a_14_next_1, a_14_next_2, r)
                    let var1 := addmod(a_14, var0, r)
                    let a_14_next_3 := calldataload(0x1944)
                    let var2 := sub(r, a_14_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_16, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_17 := calldataload(0x2224)
                    let a_15 := calldataload(0x1964)
                    let a_15_next_1 := calldataload(0x1984)
                    let a_15_next_2 := calldataload(0x19a4)
                    let var0 := mulmod(a_15_next_1, a_15_next_2, r)
                    let var1 := addmod(a_15, var0, r)
                    let a_15_next_3 := calldataload(0x19c4)
                    let var2 := sub(r, a_15_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_17, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_18 := calldataload(0x2244)
                    let a_16 := calldataload(0x19e4)
                    let a_16_next_1 := calldataload(0x1a04)
                    let a_16_next_2 := calldataload(0x1a24)
                    let var0 := mulmod(a_16_next_1, a_16_next_2, r)
                    let var1 := addmod(a_16, var0, r)
                    let a_16_next_3 := calldataload(0x1a44)
                    let var2 := sub(r, a_16_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_18, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_19 := calldataload(0x2264)
                    let a_17 := calldataload(0x1a64)
                    let a_17_next_1 := calldataload(0x1a84)
                    let a_17_next_2 := calldataload(0x1aa4)
                    let var0 := mulmod(a_17_next_1, a_17_next_2, r)
                    let var1 := addmod(a_17, var0, r)
                    let a_17_next_3 := calldataload(0x1ac4)
                    let var2 := sub(r, a_17_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_19, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_20 := calldataload(0x2284)
                    let a_18 := calldataload(0x1ae4)
                    let a_18_next_1 := calldataload(0x1b04)
                    let a_18_next_2 := calldataload(0x1b24)
                    let var0 := mulmod(a_18_next_1, a_18_next_2, r)
                    let var1 := addmod(a_18, var0, r)
                    let a_18_next_3 := calldataload(0x1b44)
                    let var2 := sub(r, a_18_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_20, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_21 := calldataload(0x22a4)
                    let a_19 := calldataload(0x1b64)
                    let a_19_next_1 := calldataload(0x1b84)
                    let a_19_next_2 := calldataload(0x1ba4)
                    let var0 := mulmod(a_19_next_1, a_19_next_2, r)
                    let var1 := addmod(a_19, var0, r)
                    let a_19_next_3 := calldataload(0x1bc4)
                    let var2 := sub(r, a_19_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_21, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_22 := calldataload(0x22c4)
                    let a_20 := calldataload(0x1be4)
                    let a_20_next_1 := calldataload(0x1c04)
                    let a_20_next_2 := calldataload(0x1c24)
                    let var0 := mulmod(a_20_next_1, a_20_next_2, r)
                    let var1 := addmod(a_20, var0, r)
                    let a_20_next_3 := calldataload(0x1c44)
                    let var2 := sub(r, a_20_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_22, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_23 := calldataload(0x22e4)
                    let a_21 := calldataload(0x1c64)
                    let a_21_next_1 := calldataload(0x1c84)
                    let a_21_next_2 := calldataload(0x1ca4)
                    let var0 := mulmod(a_21_next_1, a_21_next_2, r)
                    let var1 := addmod(a_21, var0, r)
                    let a_21_next_3 := calldataload(0x1cc4)
                    let var2 := sub(r, a_21_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_23, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_24 := calldataload(0x2304)
                    let a_22 := calldataload(0x1ce4)
                    let a_22_next_1 := calldataload(0x1d04)
                    let a_22_next_2 := calldataload(0x1d24)
                    let var0 := mulmod(a_22_next_1, a_22_next_2, r)
                    let var1 := addmod(a_22, var0, r)
                    let a_22_next_3 := calldataload(0x1d44)
                    let var2 := sub(r, a_22_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_24, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_25 := calldataload(0x2324)
                    let a_23 := calldataload(0x1d64)
                    let a_23_next_1 := calldataload(0x1d84)
                    let a_23_next_2 := calldataload(0x1da4)
                    let var0 := mulmod(a_23_next_1, a_23_next_2, r)
                    let var1 := addmod(a_23, var0, r)
                    let a_23_next_3 := calldataload(0x1dc4)
                    let var2 := sub(r, a_23_next_3)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_25, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_26 := calldataload(0x2344)
                    let a_29 := calldataload(0x1e84)
                    let c_0 := mload(0x2040)
                    let var0 := mulmod(a_29, c_0, r)
                    let a_29_next_1 := calldataload(0x1f04)
                    let var1 := addmod(var0, a_29_next_1, r)
                    let a_29_next_2 := calldataload(0x1f24)
                    let var2 := sub(r, a_29_next_2)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_26, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_27 := calldataload(0x2364)
                    let a_30 := calldataload(0x1ea4)
                    let c_0 := mload(0x2040)
                    let var0 := mulmod(a_30, c_0, r)
                    let a_30_next_1 := calldataload(0x1f44)
                    let var1 := addmod(var0, a_30_next_1, r)
                    let a_30_next_2 := calldataload(0x1f64)
                    let var2 := sub(r, a_30_next_2)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_27, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_28 := calldataload(0x2384)
                    let a_31 := calldataload(0x1ec4)
                    let c_0 := mload(0x2040)
                    let var0 := mulmod(a_31, c_0, r)
                    let a_31_next_1 := calldataload(0x1f84)
                    let var1 := addmod(var0, a_31_next_1, r)
                    let a_31_next_2 := calldataload(0x1fa4)
                    let var2 := sub(r, a_31_next_2)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_28, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let f_29 := calldataload(0x23a4)
                    let a_32 := calldataload(0x1ee4)
                    let c_0 := mload(0x2040)
                    let var0 := mulmod(a_32, c_0, r)
                    let a_32_next_1 := calldataload(0x1fc4)
                    let var1 := addmod(var0, a_32_next_1, r)
                    let a_32_next_2 := calldataload(0x1fe4)
                    let var2 := sub(r, a_32_next_2)
                    let var3 := addmod(var1, var2, r)
                    let var4 := mulmod(f_29, var3, r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), var4, r)
                }
                {
                    let l_0 := mload(L_0_MPTR)
                    let eval := addmod(l_0, sub(r, mulmod(l_0, calldataload(0x2844), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let perm_z_last := calldataload(0x2ea4)
                    let eval := mulmod(mload(L_LAST_MPTR), addmod(mulmod(perm_z_last, perm_z_last, r), sub(r, perm_z_last), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x28a4), sub(r, calldataload(0x2884)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2904), sub(r, calldataload(0x28e4)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2964), sub(r, calldataload(0x2944)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x29c4), sub(r, calldataload(0x29a4)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2a24), sub(r, calldataload(0x2a04)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2a84), sub(r, calldataload(0x2a64)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2ae4), sub(r, calldataload(0x2ac4)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2b44), sub(r, calldataload(0x2b24)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2ba4), sub(r, calldataload(0x2b84)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2c04), sub(r, calldataload(0x2be4)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2c64), sub(r, calldataload(0x2c44)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2cc4), sub(r, calldataload(0x2ca4)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2d24), sub(r, calldataload(0x2d04)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2d84), sub(r, calldataload(0x2d64)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2de4), sub(r, calldataload(0x2dc4)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2e44), sub(r, calldataload(0x2e24)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2ea4), sub(r, calldataload(0x2e84)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x2864)
                    let rhs := calldataload(0x2844)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x2004), mulmod(beta, calldataload(0x23e4), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x11e4), mulmod(beta, calldataload(0x2404), r), r), gamma, r), r)
                    mstore(0x00, mulmod(beta, mload(X_MPTR), r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x2004), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x11e4), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x28c4)
                    let rhs := calldataload(0x28a4)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1264), mulmod(beta, calldataload(0x2424), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x12e4), mulmod(beta, calldataload(0x2444), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1264), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x12e4), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x2924)
                    let rhs := calldataload(0x2904)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1364), mulmod(beta, calldataload(0x2464), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x13e4), mulmod(beta, calldataload(0x2484), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1364), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x13e4), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x2984)
                    let rhs := calldataload(0x2964)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1464), mulmod(beta, calldataload(0x24a4), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x14e4), mulmod(beta, calldataload(0x24c4), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1464), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x14e4), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x29e4)
                    let rhs := calldataload(0x29c4)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1564), mulmod(beta, calldataload(0x24e4), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x15e4), mulmod(beta, calldataload(0x2504), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1564), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x15e4), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x2a44)
                    let rhs := calldataload(0x2a24)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1664), mulmod(beta, calldataload(0x2524), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x16e4), mulmod(beta, calldataload(0x2544), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1664), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x16e4), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x2aa4)
                    let rhs := calldataload(0x2a84)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1764), mulmod(beta, calldataload(0x2564), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x17e4), mulmod(beta, calldataload(0x2584), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1764), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x17e4), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x2b04)
                    let rhs := calldataload(0x2ae4)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1864), mulmod(beta, calldataload(0x25a4), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x18e4), mulmod(beta, calldataload(0x25c4), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1864), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x18e4), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x2b64)
                    let rhs := calldataload(0x2b44)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1964), mulmod(beta, calldataload(0x25e4), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x19e4), mulmod(beta, calldataload(0x2604), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1964), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x19e4), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x2bc4)
                    let rhs := calldataload(0x2ba4)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1a64), mulmod(beta, calldataload(0x2624), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1ae4), mulmod(beta, calldataload(0x2644), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1a64), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1ae4), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x2c24)
                    let rhs := calldataload(0x2c04)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1b64), mulmod(beta, calldataload(0x2664), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1be4), mulmod(beta, calldataload(0x2684), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1b64), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1be4), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x2c84)
                    let rhs := calldataload(0x2c64)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1c64), mulmod(beta, calldataload(0x26a4), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1ce4), mulmod(beta, calldataload(0x26c4), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1c64), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1ce4), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x2ce4)
                    let rhs := calldataload(0x2cc4)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1d64), mulmod(beta, calldataload(0x26e4), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1de4), mulmod(beta, calldataload(0x2704), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1d64), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1de4), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x2d44)
                    let rhs := calldataload(0x2d24)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1e04), mulmod(beta, calldataload(0x2724), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1e24), mulmod(beta, calldataload(0x2744), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1e04), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1e24), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x2da4)
                    let rhs := calldataload(0x2d84)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1e44), mulmod(beta, calldataload(0x2764), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1e64), mulmod(beta, calldataload(0x2784), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1e44), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1e64), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x2e04)
                    let rhs := calldataload(0x2de4)
                    lhs := mulmod(lhs, addmod(addmod(mload(INSTANCE_EVAL_MPTR), mulmod(beta, calldataload(0x27a4), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1e84), mulmod(beta, calldataload(0x27c4), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(mload(INSTANCE_EVAL_MPTR), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1e84), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x2e64)
                    let rhs := calldataload(0x2e44)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1ea4), mulmod(beta, calldataload(0x27e4), r), r), gamma, r), r)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1ec4), mulmod(beta, calldataload(0x2804), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1ea4), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1ec4), mload(0x00), r), gamma, r), r)
                    mstore(0x00, mulmod(mload(0x00), delta, r))
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let gamma := mload(GAMMA_MPTR)
                    let beta := mload(BETA_MPTR)
                    let lhs := calldataload(0x2ec4)
                    let rhs := calldataload(0x2ea4)
                    lhs := mulmod(lhs, addmod(addmod(calldataload(0x1ee4), mulmod(beta, calldataload(0x2824), r), r), gamma, r), r)
                    rhs := mulmod(rhs, addmod(addmod(calldataload(0x1ee4), mload(0x00), r), gamma, r), r)
                    let left_sub_right := addmod(lhs, sub(r, rhs), r)
                    let eval := addmod(left_sub_right, sub(r, mulmod(left_sub_right, addmod(mload(L_LAST_MPTR), mload(L_BLIND_MPTR), r), r)), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let l_0 := mload(L_0_MPTR)
                    let eval := addmod(l_0, mulmod(l_0, sub(r, calldataload(0x2ee4)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let l_last := mload(L_LAST_MPTR)
                    let eval := mulmod(l_last, addmod(mulmod(calldataload(0x2ee4), calldataload(0x2ee4), r), sub(r, calldataload(0x2ee4)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let theta := mload(THETA_MPTR)
                    let input
                    {
                        let a_24 := calldataload(0x1de4)
                        input := a_24
                    }
                    let table
                    {
                        let f_0 := calldataload(0x2024)
                        table := f_0
                    }
                    let beta := mload(BETA_MPTR)
                    let gamma := mload(GAMMA_MPTR)
                    let lhs := mulmod(calldataload(0x2f04), mulmod(addmod(calldataload(0x2f24), beta, r), addmod(calldataload(0x2f64), gamma, r), r), r)
                    let rhs := mulmod(calldataload(0x2ee4), mulmod(addmod(input, beta, r), addmod(table, gamma, r), r), r)
                    let eval := mulmod(addmod(1, sub(r, addmod(mload(L_BLIND_MPTR), mload(L_LAST_MPTR), r)), r), addmod(lhs, sub(r, rhs), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2f24), sub(r, calldataload(0x2f64)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(addmod(1, sub(r, addmod(mload(L_BLIND_MPTR), mload(L_LAST_MPTR), r)), r), mulmod(addmod(calldataload(0x2f24), sub(r, calldataload(0x2f64)), r), addmod(calldataload(0x2f24), sub(r, calldataload(0x2f44)), r), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let l_0 := mload(L_0_MPTR)
                    let eval := addmod(l_0, mulmod(l_0, sub(r, calldataload(0x2f84)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let l_last := mload(L_LAST_MPTR)
                    let eval := mulmod(l_last, addmod(mulmod(calldataload(0x2f84), calldataload(0x2f84), r), sub(r, calldataload(0x2f84)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let theta := mload(THETA_MPTR)
                    let input
                    {
                        let a_25 := calldataload(0x1e04)
                        input := a_25
                    }
                    let table
                    {
                        let f_0 := calldataload(0x2024)
                        table := f_0
                    }
                    let beta := mload(BETA_MPTR)
                    let gamma := mload(GAMMA_MPTR)
                    let lhs := mulmod(calldataload(0x2fa4), mulmod(addmod(calldataload(0x2fc4), beta, r), addmod(calldataload(0x3004), gamma, r), r), r)
                    let rhs := mulmod(calldataload(0x2f84), mulmod(addmod(input, beta, r), addmod(table, gamma, r), r), r)
                    let eval := mulmod(addmod(1, sub(r, addmod(mload(L_BLIND_MPTR), mload(L_LAST_MPTR), r)), r), addmod(lhs, sub(r, rhs), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x2fc4), sub(r, calldataload(0x3004)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(addmod(1, sub(r, addmod(mload(L_BLIND_MPTR), mload(L_LAST_MPTR), r)), r), mulmod(addmod(calldataload(0x2fc4), sub(r, calldataload(0x3004)), r), addmod(calldataload(0x2fc4), sub(r, calldataload(0x2fe4)), r), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let l_0 := mload(L_0_MPTR)
                    let eval := addmod(l_0, mulmod(l_0, sub(r, calldataload(0x3024)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let l_last := mload(L_LAST_MPTR)
                    let eval := mulmod(l_last, addmod(mulmod(calldataload(0x3024), calldataload(0x3024), r), sub(r, calldataload(0x3024)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let theta := mload(THETA_MPTR)
                    let input
                    {
                        let a_26 := calldataload(0x1e24)
                        input := a_26
                    }
                    let table
                    {
                        let f_0 := calldataload(0x2024)
                        table := f_0
                    }
                    let beta := mload(BETA_MPTR)
                    let gamma := mload(GAMMA_MPTR)
                    let lhs := mulmod(calldataload(0x3044), mulmod(addmod(calldataload(0x3064), beta, r), addmod(calldataload(0x30a4), gamma, r), r), r)
                    let rhs := mulmod(calldataload(0x3024), mulmod(addmod(input, beta, r), addmod(table, gamma, r), r), r)
                    let eval := mulmod(addmod(1, sub(r, addmod(mload(L_BLIND_MPTR), mload(L_LAST_MPTR), r)), r), addmod(lhs, sub(r, rhs), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x3064), sub(r, calldataload(0x30a4)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(addmod(1, sub(r, addmod(mload(L_BLIND_MPTR), mload(L_LAST_MPTR), r)), r), mulmod(addmod(calldataload(0x3064), sub(r, calldataload(0x30a4)), r), addmod(calldataload(0x3064), sub(r, calldataload(0x3084)), r), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let l_0 := mload(L_0_MPTR)
                    let eval := addmod(l_0, mulmod(l_0, sub(r, calldataload(0x30c4)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let l_last := mload(L_LAST_MPTR)
                    let eval := mulmod(l_last, addmod(mulmod(calldataload(0x30c4), calldataload(0x30c4), r), sub(r, calldataload(0x30c4)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let theta := mload(THETA_MPTR)
                    let input
                    {
                        let a_27 := calldataload(0x1e44)
                        input := a_27
                    }
                    let table
                    {
                        let f_0 := calldataload(0x2024)
                        table := f_0
                    }
                    let beta := mload(BETA_MPTR)
                    let gamma := mload(GAMMA_MPTR)
                    let lhs := mulmod(calldataload(0x30e4), mulmod(addmod(calldataload(0x3104), beta, r), addmod(calldataload(0x3144), gamma, r), r), r)
                    let rhs := mulmod(calldataload(0x30c4), mulmod(addmod(input, beta, r), addmod(table, gamma, r), r), r)
                    let eval := mulmod(addmod(1, sub(r, addmod(mload(L_BLIND_MPTR), mload(L_LAST_MPTR), r)), r), addmod(lhs, sub(r, rhs), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x3104), sub(r, calldataload(0x3144)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(addmod(1, sub(r, addmod(mload(L_BLIND_MPTR), mload(L_LAST_MPTR), r)), r), mulmod(addmod(calldataload(0x3104), sub(r, calldataload(0x3144)), r), addmod(calldataload(0x3104), sub(r, calldataload(0x3124)), r), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let l_0 := mload(L_0_MPTR)
                    let eval := addmod(l_0, mulmod(l_0, sub(r, calldataload(0x3164)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let l_last := mload(L_LAST_MPTR)
                    let eval := mulmod(l_last, addmod(mulmod(calldataload(0x3164), calldataload(0x3164), r), sub(r, calldataload(0x3164)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let theta := mload(THETA_MPTR)
                    let input
                    {
                        let a_28 := calldataload(0x1e64)
                        input := a_28
                    }
                    let table
                    {
                        let f_0 := calldataload(0x2024)
                        table := f_0
                    }
                    let beta := mload(BETA_MPTR)
                    let gamma := mload(GAMMA_MPTR)
                    let lhs := mulmod(calldataload(0x3184), mulmod(addmod(calldataload(0x31a4), beta, r), addmod(calldataload(0x31e4), gamma, r), r), r)
                    let rhs := mulmod(calldataload(0x3164), mulmod(addmod(input, beta, r), addmod(table, gamma, r), r), r)
                    let eval := mulmod(addmod(1, sub(r, addmod(mload(L_BLIND_MPTR), mload(L_LAST_MPTR), r)), r), addmod(lhs, sub(r, rhs), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(mload(L_0_MPTR), addmod(calldataload(0x31a4), sub(r, calldataload(0x31e4)), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }
                {
                    let eval := mulmod(addmod(1, sub(r, addmod(mload(L_BLIND_MPTR), mload(L_LAST_MPTR), r)), r), mulmod(addmod(calldataload(0x31a4), sub(r, calldataload(0x31e4)), r), addmod(calldataload(0x31a4), sub(r, calldataload(0x31c4)), r), r), r)
                    quotient_eval_numer := addmod(mulmod(quotient_eval_numer, y, r), eval, r)
                }

                pop(y)
                pop(delta)

                let quotient_eval := mulmod(quotient_eval_numer, mload(X_N_MINUS_1_INV_MPTR), r)
                mstore(QUOTIENT_EVAL_MPTR, quotient_eval)
            }

            // Compute quotient commitment
            {
                mstore(0x00, calldataload(LAST_QUOTIENT_X_CPTR))
                mstore(0x20, calldataload(add(LAST_QUOTIENT_X_CPTR, 0x20)))
                let x_n := mload(X_N_MPTR)
                for
                    {
                        let cptr := sub(LAST_QUOTIENT_X_CPTR, 0x40)
                        let cptr_end := sub(FIRST_QUOTIENT_X_CPTR, 0x40)
                    }
                    lt(cptr_end, cptr)
                    {}
                {
                    success := ec_mul_acc(success, x_n)
                    success := ec_add_acc(success, calldataload(cptr), calldataload(add(cptr, 0x20)))
                    cptr := sub(cptr, 0x40)
                }
                mstore(QUOTIENT_X_MPTR, mload(0x00))
                mstore(QUOTIENT_Y_MPTR, mload(0x20))
            }

            // Compute pairing lhs and rhs
            {
                {
                    let x := mload(X_MPTR)
                    let omega := mload(OMEGA_MPTR)
                    let omega_inv := mload(OMEGA_INV_MPTR)
                    let x_pow_of_omega := mulmod(x, omega, r)
                    mstore(0x0520, x_pow_of_omega)
                    x_pow_of_omega := mulmod(x_pow_of_omega, omega, r)
                    mstore(0x0540, x_pow_of_omega)
                    x_pow_of_omega := mulmod(x_pow_of_omega, omega, r)
                    mstore(0x0560, x_pow_of_omega)
                    mstore(0x0500, x)
                    x_pow_of_omega := mulmod(x, omega_inv, r)
                    mstore(0x04e0, x_pow_of_omega)
                    x_pow_of_omega := mulmod(x_pow_of_omega, omega_inv, r)
                    x_pow_of_omega := mulmod(x_pow_of_omega, omega_inv, r)
                    x_pow_of_omega := mulmod(x_pow_of_omega, omega_inv, r)
                    x_pow_of_omega := mulmod(x_pow_of_omega, omega_inv, r)
                    x_pow_of_omega := mulmod(x_pow_of_omega, omega_inv, r)
                    x_pow_of_omega := mulmod(x_pow_of_omega, omega_inv, r)
                    mstore(0x04c0, x_pow_of_omega)
                }
                {
                    let mu := mload(MU_MPTR)
                    for
                        {
                            let mptr := 0x0580
                            let mptr_end := 0x0640
                            let point_mptr := 0x04c0
                        }
                        lt(mptr, mptr_end)
                        {
                            mptr := add(mptr, 0x20)
                            point_mptr := add(point_mptr, 0x20)
                        }
                    {
                        mstore(mptr, addmod(mu, sub(r, mload(point_mptr)), r))
                    }
                    let s
                    s := mload(0x05c0)
                    s := mulmod(s, mload(0x05e0), r)
                    s := mulmod(s, mload(0x0600), r)
                    s := mulmod(s, mload(0x0620), r)
                    mstore(0x0640, s)
                    let diff
                    diff := mload(0x0580)
                    diff := mulmod(diff, mload(0x05a0), r)
                    mstore(0x0660, diff)
                    mstore(0x00, diff)
                    diff := mload(0x0580)
                    diff := mulmod(diff, mload(0x05a0), r)
                    diff := mulmod(diff, mload(0x05e0), r)
                    diff := mulmod(diff, mload(0x0600), r)
                    diff := mulmod(diff, mload(0x0620), r)
                    mstore(0x0680, diff)
                    diff := mload(0x0580)
                    diff := mulmod(diff, mload(0x05a0), r)
                    diff := mulmod(diff, mload(0x0620), r)
                    mstore(0x06a0, diff)
                    diff := mload(0x05a0)
                    diff := mulmod(diff, mload(0x0600), r)
                    diff := mulmod(diff, mload(0x0620), r)
                    mstore(0x06c0, diff)
                    diff := mload(0x0580)
                    diff := mulmod(diff, mload(0x05a0), r)
                    diff := mulmod(diff, mload(0x0600), r)
                    diff := mulmod(diff, mload(0x0620), r)
                    mstore(0x06e0, diff)
                    diff := mload(0x0580)
                    diff := mulmod(diff, mload(0x05e0), r)
                    diff := mulmod(diff, mload(0x0600), r)
                    diff := mulmod(diff, mload(0x0620), r)
                    mstore(0x0700, diff)
                }
                {
                    let point_2 := mload(0x0500)
                    let point_3 := mload(0x0520)
                    let point_4 := mload(0x0540)
                    let point_5 := mload(0x0560)
                    let coeff
                    coeff := addmod(point_2, sub(r, point_3), r)
                    coeff := mulmod(coeff, addmod(point_2, sub(r, point_4), r), r)
                    coeff := mulmod(coeff, addmod(point_2, sub(r, point_5), r), r)
                    coeff := mulmod(coeff, mload(0x05c0), r)
                    mstore(0x20, coeff)
                    coeff := addmod(point_3, sub(r, point_2), r)
                    coeff := mulmod(coeff, addmod(point_3, sub(r, point_4), r), r)
                    coeff := mulmod(coeff, addmod(point_3, sub(r, point_5), r), r)
                    coeff := mulmod(coeff, mload(0x05e0), r)
                    mstore(0x40, coeff)
                    coeff := addmod(point_4, sub(r, point_2), r)
                    coeff := mulmod(coeff, addmod(point_4, sub(r, point_3), r), r)
                    coeff := mulmod(coeff, addmod(point_4, sub(r, point_5), r), r)
                    coeff := mulmod(coeff, mload(0x0600), r)
                    mstore(0x60, coeff)
                    coeff := addmod(point_5, sub(r, point_2), r)
                    coeff := mulmod(coeff, addmod(point_5, sub(r, point_3), r), r)
                    coeff := mulmod(coeff, addmod(point_5, sub(r, point_4), r), r)
                    coeff := mulmod(coeff, mload(0x0620), r)
                    mstore(0x80, coeff)
                }
                {
                    let point_2 := mload(0x0500)
                    let coeff
                    coeff := 1
                    coeff := mulmod(coeff, mload(0x05c0), r)
                    mstore(0xa0, coeff)
                }
                {
                    let point_2 := mload(0x0500)
                    let point_3 := mload(0x0520)
                    let point_4 := mload(0x0540)
                    let coeff
                    coeff := addmod(point_2, sub(r, point_3), r)
                    coeff := mulmod(coeff, addmod(point_2, sub(r, point_4), r), r)
                    coeff := mulmod(coeff, mload(0x05c0), r)
                    mstore(0xc0, coeff)
                    coeff := addmod(point_3, sub(r, point_2), r)
                    coeff := mulmod(coeff, addmod(point_3, sub(r, point_4), r), r)
                    coeff := mulmod(coeff, mload(0x05e0), r)
                    mstore(0xe0, coeff)
                    coeff := addmod(point_4, sub(r, point_2), r)
                    coeff := mulmod(coeff, addmod(point_4, sub(r, point_3), r), r)
                    coeff := mulmod(coeff, mload(0x0600), r)
                    mstore(0x0100, coeff)
                }
                {
                    let point_0 := mload(0x04c0)
                    let point_2 := mload(0x0500)
                    let point_3 := mload(0x0520)
                    let coeff
                    coeff := addmod(point_0, sub(r, point_2), r)
                    coeff := mulmod(coeff, addmod(point_0, sub(r, point_3), r), r)
                    coeff := mulmod(coeff, mload(0x0580), r)
                    mstore(0x0120, coeff)
                    coeff := addmod(point_2, sub(r, point_0), r)
                    coeff := mulmod(coeff, addmod(point_2, sub(r, point_3), r), r)
                    coeff := mulmod(coeff, mload(0x05c0), r)
                    mstore(0x0140, coeff)
                    coeff := addmod(point_3, sub(r, point_0), r)
                    coeff := mulmod(coeff, addmod(point_3, sub(r, point_2), r), r)
                    coeff := mulmod(coeff, mload(0x05e0), r)
                    mstore(0x0160, coeff)
                }
                {
                    let point_2 := mload(0x0500)
                    let point_3 := mload(0x0520)
                    let coeff
                    coeff := addmod(point_2, sub(r, point_3), r)
                    coeff := mulmod(coeff, mload(0x05c0), r)
                    mstore(0x0180, coeff)
                    coeff := addmod(point_3, sub(r, point_2), r)
                    coeff := mulmod(coeff, mload(0x05e0), r)
                    mstore(0x01a0, coeff)
                }
                {
                    let point_1 := mload(0x04e0)
                    let point_2 := mload(0x0500)
                    let coeff
                    coeff := addmod(point_1, sub(r, point_2), r)
                    coeff := mulmod(coeff, mload(0x05a0), r)
                    mstore(0x01c0, coeff)
                    coeff := addmod(point_2, sub(r, point_1), r)
                    coeff := mulmod(coeff, mload(0x05c0), r)
                    mstore(0x01e0, coeff)
                }
                {
                    success := batch_invert(success, 0, 0x0200, r)
                    let diff_0_inv := mload(0x00)
                    mstore(0x0660, diff_0_inv)
                    for
                        {
                            let mptr := 0x0680
                            let mptr_end := 0x0720
                        }
                        lt(mptr, mptr_end)
                        { mptr := add(mptr, 0x20) }
                    {
                        mstore(mptr, mulmod(mload(mptr), diff_0_inv, r))
                    }
                }
                {
                    let zeta := mload(ZETA_MPTR)
                    let r_eval
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x1d64), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1d84), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x1da4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x1dc4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x1ce4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1d04), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x1d24), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x1d44), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x1c64), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1c84), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x1ca4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x1cc4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x1be4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1c04), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x1c24), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x1c44), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x1b64), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1b84), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x1ba4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x1bc4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x1ae4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1b04), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x1b24), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x1b44), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x1a64), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1a84), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x1aa4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x1ac4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x19e4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1a04), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x1a24), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x1a44), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x1964), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1984), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x19a4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x19c4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x18e4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1904), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x1924), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x1944), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x1864), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1884), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x18a4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x18c4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x17e4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1804), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x1824), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x1844), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x1764), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1784), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x17a4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x17c4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x16e4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1704), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x1724), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x1744), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x1664), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1684), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x16a4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x16c4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x15e4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1604), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x1624), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x1644), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x1564), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1584), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x15a4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x15c4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x14e4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1504), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x1524), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x1544), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x1464), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1484), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x14a4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x14c4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x13e4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1404), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x1424), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x1444), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x1364), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1384), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x13a4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x13c4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x12e4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1304), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x1324), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x1344), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x1264), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1284), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x12a4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x12c4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x20), calldataload(0x11e4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x40), calldataload(0x1204), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x60), calldataload(0x1224), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x80), calldataload(0x1244), r), r)
                    mstore(0x0720, r_eval)
                }
                {
                    let coeff := mload(0xa0)
                    let zeta := mload(ZETA_MPTR)
                    let r_eval
                    r_eval := mulmod(coeff, calldataload(0x23c4), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(coeff, mload(QUOTIENT_EVAL_MPTR), r), r)
                    for
                        {
                            let cptr := 0x2824
                            let cptr_end := 0x23c4
                        }
                        lt(cptr_end, cptr)
                        { cptr := sub(cptr, 0x20) }
                    {
                        r_eval := addmod(mulmod(r_eval, zeta, r), mulmod(coeff, calldataload(cptr), r), r)
                    }
                    for
                        {
                            let cptr := 0x23a4
                            let cptr_end := 0x1fe4
                        }
                        lt(cptr_end, cptr)
                        { cptr := sub(cptr, 0x20) }
                    {
                        r_eval := addmod(mulmod(r_eval, zeta, r), mulmod(coeff, calldataload(cptr), r), r)
                    }
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(coeff, calldataload(0x31e4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(coeff, calldataload(0x3144), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(coeff, calldataload(0x30a4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(coeff, calldataload(0x3004), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(coeff, calldataload(0x2f64), r), r)
                    for
                        {
                            let cptr := 0x1e64
                            let cptr_end := 0x1dc4
                        }
                        lt(cptr_end, cptr)
                        { cptr := sub(cptr, 0x20) }
                    {
                        r_eval := addmod(mulmod(r_eval, zeta, r), mulmod(coeff, calldataload(cptr), r), r)
                    }
                    r_eval := mulmod(r_eval, mload(0x0680), r)
                    mstore(0x0740, r_eval)
                }
                {
                    let zeta := mload(ZETA_MPTR)
                    let r_eval
                    r_eval := addmod(r_eval, mulmod(mload(0xc0), calldataload(0x1ee4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0xe0), calldataload(0x1fc4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0100), calldataload(0x1fe4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0xc0), calldataload(0x1ec4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0xe0), calldataload(0x1f84), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0100), calldataload(0x1fa4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0xc0), calldataload(0x1ea4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0xe0), calldataload(0x1f44), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0100), calldataload(0x1f64), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0xc0), calldataload(0x1e84), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0xe0), calldataload(0x1f04), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0100), calldataload(0x1f24), r), r)
                    r_eval := mulmod(r_eval, mload(0x06a0), r)
                    mstore(0x0760, r_eval)
                }
                {
                    let zeta := mload(ZETA_MPTR)
                    let r_eval
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x2e84), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x2e44), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x2e64), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x2e24), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x2de4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x2e04), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x2dc4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x2d84), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x2da4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x2d64), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x2d24), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x2d44), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x2d04), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x2cc4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x2ce4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x2ca4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x2c64), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x2c84), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x2c44), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x2c04), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x2c24), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x2be4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x2ba4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x2bc4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x2b84), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x2b44), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x2b64), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x2b24), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x2ae4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x2b04), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x2ac4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x2a84), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x2aa4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x2a64), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x2a24), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x2a44), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x2a04), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x29c4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x29e4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x29a4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x2964), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x2984), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x2944), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x2904), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x2924), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x28e4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x28a4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x28c4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0120), calldataload(0x2884), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0140), calldataload(0x2844), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0160), calldataload(0x2864), r), r)
                    r_eval := mulmod(r_eval, mload(0x06c0), r)
                    mstore(0x0780, r_eval)
                }
                {
                    let zeta := mload(ZETA_MPTR)
                    let r_eval
                    r_eval := addmod(r_eval, mulmod(mload(0x0180), calldataload(0x3164), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x01a0), calldataload(0x3184), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0180), calldataload(0x30c4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x01a0), calldataload(0x30e4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0180), calldataload(0x3024), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x01a0), calldataload(0x3044), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0180), calldataload(0x2f84), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x01a0), calldataload(0x2fa4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0180), calldataload(0x2ee4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x01a0), calldataload(0x2f04), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x0180), calldataload(0x2ea4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x01a0), calldataload(0x2ec4), r), r)
                    r_eval := mulmod(r_eval, mload(0x06e0), r)
                    mstore(0x07a0, r_eval)
                }
                {
                    let zeta := mload(ZETA_MPTR)
                    let r_eval
                    r_eval := addmod(r_eval, mulmod(mload(0x01c0), calldataload(0x31c4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x01e0), calldataload(0x31a4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x01c0), calldataload(0x3124), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x01e0), calldataload(0x3104), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x01c0), calldataload(0x3084), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x01e0), calldataload(0x3064), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x01c0), calldataload(0x2fe4), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x01e0), calldataload(0x2fc4), r), r)
                    r_eval := mulmod(r_eval, zeta, r)
                    r_eval := addmod(r_eval, mulmod(mload(0x01c0), calldataload(0x2f44), r), r)
                    r_eval := addmod(r_eval, mulmod(mload(0x01e0), calldataload(0x2f24), r), r)
                    r_eval := mulmod(r_eval, mload(0x0700), r)
                    mstore(0x07c0, r_eval)
                }
                {
                    let sum := mload(0x20)
                    sum := addmod(sum, mload(0x40), r)
                    sum := addmod(sum, mload(0x60), r)
                    sum := addmod(sum, mload(0x80), r)
                    mstore(0x07e0, sum)
                }
                {
                    let sum := mload(0xa0)
                    mstore(0x0800, sum)
                }
                {
                    let sum := mload(0xc0)
                    sum := addmod(sum, mload(0xe0), r)
                    sum := addmod(sum, mload(0x0100), r)
                    mstore(0x0820, sum)
                }
                {
                    let sum := mload(0x0120)
                    sum := addmod(sum, mload(0x0140), r)
                    sum := addmod(sum, mload(0x0160), r)
                    mstore(0x0840, sum)
                }
                {
                    let sum := mload(0x0180)
                    sum := addmod(sum, mload(0x01a0), r)
                    mstore(0x0860, sum)
                }
                {
                    let sum := mload(0x01c0)
                    sum := addmod(sum, mload(0x01e0), r)
                    mstore(0x0880, sum)
                }
                {
                    for
                        {
                            let mptr := 0x00
                            let mptr_end := 0xc0
                            let sum_mptr := 0x07e0
                        }
                        lt(mptr, mptr_end)
                        {
                            mptr := add(mptr, 0x20)
                            sum_mptr := add(sum_mptr, 0x20)
                        }
                    {
                        mstore(mptr, mload(sum_mptr))
                    }
                    success := batch_invert(success, 0, 0xc0, r)
                    let r_eval := mulmod(mload(0xa0), mload(0x07c0), r)
                    for
                        {
                            let sum_inv_mptr := 0x80
                            let sum_inv_mptr_end := 0xc0
                            let r_eval_mptr := 0x07a0
                        }
                        lt(sum_inv_mptr, sum_inv_mptr_end)
                        {
                            sum_inv_mptr := sub(sum_inv_mptr, 0x20)
                            r_eval_mptr := sub(r_eval_mptr, 0x20)
                        }
                    {
                        r_eval := mulmod(r_eval, mload(NU_MPTR), r)
                        r_eval := addmod(r_eval, mulmod(mload(sum_inv_mptr), mload(r_eval_mptr), r), r)
                    }
                    mstore(G1_SCALAR_MPTR, sub(r, r_eval))
                }
                {
                    let zeta := mload(ZETA_MPTR)
                    let nu := mload(NU_MPTR)
                    mstore(0x00, calldataload(0x0624))
                    mstore(0x20, calldataload(0x0644))
                    for
                        {
                            let ptr := 0x05e4
                            let ptr_end := 0x24
                        }
                        lt(ptr_end, ptr)
                        { ptr := sub(ptr, 0x40) }
                    {
                        success := ec_mul_acc(success, zeta)
                        success := ec_add_acc(success, calldataload(ptr), calldataload(add(ptr, 0x20)))
                    }
                    mstore(0x80, calldataload(0x10e4))
                    mstore(0xa0, calldataload(0x1104))
                    success := ec_mul_tmp(success, zeta)
                    success := ec_add_tmp(success, mload(QUOTIENT_X_MPTR), mload(QUOTIENT_Y_MPTR))
                    for
                        {
                            let ptr := 0x2000
                            let ptr_end := 0x1040
                        }
                        lt(ptr_end, ptr)
                        { ptr := sub(ptr, 0x40) }
                    {
                        success := ec_mul_tmp(success, zeta)
                        success := ec_add_tmp(success, mload(ptr), mload(add(ptr, 0x20)))
                    }
                    success := ec_mul_tmp(success, zeta)
                    success := ec_add_tmp(success, mload(0x1000), mload(0x1020))
                    success := ec_mul_tmp(success, zeta)
                    success := ec_add_tmp(success, mload(0x1040), mload(0x1060))
                    success := ec_mul_tmp(success, zeta)
                    success := ec_add_tmp(success, calldataload(0x0ae4), calldataload(0x0b04))
                    success := ec_mul_tmp(success, zeta)
                    success := ec_add_tmp(success, calldataload(0x0a64), calldataload(0x0a84))
                    success := ec_mul_tmp(success, zeta)
                    success := ec_add_tmp(success, calldataload(0x09e4), calldataload(0x0a04))
                    success := ec_mul_tmp(success, zeta)
                    success := ec_add_tmp(success, calldataload(0x0964), calldataload(0x0984))
                    success := ec_mul_tmp(success, zeta)
                    success := ec_add_tmp(success, calldataload(0x08e4), calldataload(0x0904))
                    for
                        {
                            let ptr := 0x0764
                            let ptr_end := 0x0624
                        }
                        lt(ptr_end, ptr)
                        { ptr := sub(ptr, 0x40) }
                    {
                        success := ec_mul_tmp(success, zeta)
                        success := ec_add_tmp(success, calldataload(ptr), calldataload(add(ptr, 0x20)))
                    }
                    success := ec_mul_tmp(success, mulmod(nu, mload(0x0680), r))
                    success := ec_add_acc(success, mload(0x80), mload(0xa0))
                    nu := mulmod(nu, mload(NU_MPTR), r)
                    mstore(0x80, calldataload(0x0864))
                    mstore(0xa0, calldataload(0x0884))
                    for
                        {
                            let ptr := 0x0824
                            let ptr_end := 0x0764
                        }
                        lt(ptr_end, ptr)
                        { ptr := sub(ptr, 0x40) }
                    {
                        success := ec_mul_tmp(success, zeta)
                        success := ec_add_tmp(success, calldataload(ptr), calldataload(add(ptr, 0x20)))
                    }
                    success := ec_mul_tmp(success, mulmod(nu, mload(0x06a0), r))
                    success := ec_add_acc(success, mload(0x80), mload(0xa0))
                    nu := mulmod(nu, mload(NU_MPTR), r)
                    mstore(0x80, calldataload(0x0f24))
                    mstore(0xa0, calldataload(0x0f44))
                    for
                        {
                            let ptr := 0x0ee4
                            let ptr_end := 0x0ae4
                        }
                        lt(ptr_end, ptr)
                        { ptr := sub(ptr, 0x40) }
                    {
                        success := ec_mul_tmp(success, zeta)
                        success := ec_add_tmp(success, calldataload(ptr), calldataload(add(ptr, 0x20)))
                    }
                    success := ec_mul_tmp(success, mulmod(nu, mload(0x06c0), r))
                    success := ec_add_acc(success, mload(0x80), mload(0xa0))
                    nu := mulmod(nu, mload(NU_MPTR), r)
                    mstore(0x80, calldataload(0x10a4))
                    mstore(0xa0, calldataload(0x10c4))
                    for
                        {
                            let ptr := 0x1064
                            let ptr_end := 0x0f24
                        }
                        lt(ptr_end, ptr)
                        { ptr := sub(ptr, 0x40) }
                    {
                        success := ec_mul_tmp(success, zeta)
                        success := ec_add_tmp(success, calldataload(ptr), calldataload(add(ptr, 0x20)))
                    }
                    success := ec_mul_tmp(success, mulmod(nu, mload(0x06e0), r))
                    success := ec_add_acc(success, mload(0x80), mload(0xa0))
                    nu := mulmod(nu, mload(NU_MPTR), r)
                    mstore(0x80, calldataload(0x0aa4))
                    mstore(0xa0, calldataload(0x0ac4))
                    success := ec_mul_tmp(success, zeta)
                    success := ec_add_tmp(success, calldataload(0x0a24), calldataload(0x0a44))
                    success := ec_mul_tmp(success, zeta)
                    success := ec_add_tmp(success, calldataload(0x09a4), calldataload(0x09c4))
                    success := ec_mul_tmp(success, zeta)
                    success := ec_add_tmp(success, calldataload(0x0924), calldataload(0x0944))
                    success := ec_mul_tmp(success, zeta)
                    success := ec_add_tmp(success, calldataload(0x08a4), calldataload(0x08c4))
                    success := ec_mul_tmp(success, mulmod(nu, mload(0x0700), r))
                    success := ec_add_acc(success, mload(0x80), mload(0xa0))
                    mstore(0x80, mload(G1_X_MPTR))
                    mstore(0xa0, mload(G1_Y_MPTR))
                    success := ec_mul_tmp(success, mload(G1_SCALAR_MPTR))
                    success := ec_add_acc(success, mload(0x80), mload(0xa0))
                    mstore(0x80, calldataload(0x3204))
                    mstore(0xa0, calldataload(0x3224))
                    success := ec_mul_tmp(success, sub(r, mload(0x0640)))
                    success := ec_add_acc(success, mload(0x80), mload(0xa0))
                    mstore(0x80, calldataload(0x3244))
                    mstore(0xa0, calldataload(0x3264))
                    success := ec_mul_tmp(success, mload(MU_MPTR))
                    success := ec_add_acc(success, mload(0x80), mload(0xa0))
                    mstore(PAIRING_LHS_X_MPTR, mload(0x00))
                    mstore(PAIRING_LHS_Y_MPTR, mload(0x20))
                    mstore(PAIRING_RHS_X_MPTR, calldataload(0x3244))
                    mstore(PAIRING_RHS_Y_MPTR, calldataload(0x3264))
                }
            }

            // Random linear combine with accumulator
            if mload(HAS_ACCUMULATOR_MPTR) {
                mstore(0x00, mload(ACC_LHS_X_MPTR))
                mstore(0x20, mload(ACC_LHS_Y_MPTR))
                mstore(0x40, mload(ACC_RHS_X_MPTR))
                mstore(0x60, mload(ACC_RHS_Y_MPTR))
                mstore(0x80, mload(PAIRING_LHS_X_MPTR))
                mstore(0xa0, mload(PAIRING_LHS_Y_MPTR))
                mstore(0xc0, mload(PAIRING_RHS_X_MPTR))
                mstore(0xe0, mload(PAIRING_RHS_Y_MPTR))
                let challenge := mod(keccak256(0x00, 0x100), r)

                // [pairing_lhs] += challenge * [acc_lhs]
                success := ec_mul_acc(success, challenge)
                success := ec_add_acc(success, mload(PAIRING_LHS_X_MPTR), mload(PAIRING_LHS_Y_MPTR))
                mstore(PAIRING_LHS_X_MPTR, mload(0x00))
                mstore(PAIRING_LHS_Y_MPTR, mload(0x20))

                // [pairing_rhs] += challenge * [acc_rhs]
                mstore(0x00, mload(ACC_RHS_X_MPTR))
                mstore(0x20, mload(ACC_RHS_Y_MPTR))
                success := ec_mul_acc(success, challenge)
                success := ec_add_acc(success, mload(PAIRING_RHS_X_MPTR), mload(PAIRING_RHS_Y_MPTR))
                mstore(PAIRING_RHS_X_MPTR, mload(0x00))
                mstore(PAIRING_RHS_Y_MPTR, mload(0x20))
            }

            // Perform pairing
            success := ec_pairing(
                success,
                mload(PAIRING_LHS_X_MPTR),
                mload(PAIRING_LHS_Y_MPTR),
                mload(PAIRING_RHS_X_MPTR),
                mload(PAIRING_RHS_Y_MPTR)
            )

            // Revert if anything fails
            if iszero(success) {
                revert(0x00, 0x00)
            }

            // Return 1 as result if everything succeeds
            mstore(0x00, 1)
            return(0x00, 0x20)
        }
    }
}