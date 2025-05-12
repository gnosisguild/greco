import { UltraHonkBackend } from "@aztec/bb.js";
import fs from "fs";
import circuit from "../target/testing.json" assert { type: "json" };
// @ts-ignore
import { Noir } from "@noir-lang/noir_js";

// Add type for the circuit
interface NoirCircuit {
    bytecode: string;
}

const noirCircuit = circuit as unknown as NoirCircuit;

// Helper function to create a polynomial with proper structure
function createPolynomial(size: number): { coefficients: number[] } {
    return {
        coefficients: Array(size).fill(0)
    };
}

// Helper function to create an array of polynomials
function createPolynomialArray(size: number, polySize: number): { coefficients: number[] }[] {
    return Array(size).fill(null).map(() => createPolynomial(polySize));
}

(async () => {
    try {
        console.log("Initializing Noir and backend...");
        const noir = new Noir(circuit as any);
        const honk = new UltraHonkBackend(noirCircuit.bytecode, { threads: 4 });

        console.log("Creating circuit inputs...");
        // Create default inputs for the BFV encryption circuit
        const inputs = {
            pk0is: createPolynomialArray(2, 1024), // L=2, N=1024
            pk1is: createPolynomialArray(2, 1024),
            ct0is: createPolynomialArray(2, 1024),
            ct1is: createPolynomialArray(2, 1024),
            u: createPolynomial(1024),
            e0: createPolynomial(1024),
            e1: createPolynomial(1024),
            k1: createPolynomial(1024),
            r1is: createPolynomialArray(2, 2047), // (2*N)-1 = 2047
            r2is: createPolynomialArray(2, 1023), // N-1 = 1023
            p1is: createPolynomialArray(2, 2047),
            p2is: createPolynomialArray(2, 1023)
        };

        console.log("Generating witness...");
        const { witness } = await noir.execute(inputs);

        console.log("Generating proof (this may take a while for the large circuit)...");
        const { proof, publicInputs } = await honk.generateProof(witness);

        console.log("Saving proof and public inputs...");
        fs.writeFileSync("../target/proof", proof);
        fs.writeFileSync("../target/public-inputs", JSON.stringify(publicInputs));

        console.log("Proof generated successfully");
        process.exit(0);
    } catch (error) {
        console.error("Error details:", error);
        if (error instanceof Error) {
            console.error("Error stack:", error.stack);
        }
        process.exit(1);
    }
})();
