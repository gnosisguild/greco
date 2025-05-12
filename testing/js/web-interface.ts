import { UltraHonkBackend } from "@aztec/bb.js";
import { Noir } from "@noir-lang/noir_js";
import circuit from "../target/testing.json";

let currentProof: Uint8Array | null = null;
let currentPublicInputs: any = null;

const generateBtn = document.getElementById('generateBtn') as HTMLButtonElement;
const verifyBtn = document.getElementById('verifyBtn') as HTMLButtonElement;
const statusDiv = document.getElementById('status') as HTMLDivElement;

function updateStatus(message: string, isError: boolean = false) {
    console.log('Status update:', message);
    statusDiv.textContent = message;
    statusDiv.className = isError ? 'error' : 'success';
    statusDiv.style.display = 'block';
}

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

interface NoirCircuit {
    bytecode: string;
}

const noirCircuit = circuit as unknown as NoirCircuit;

async function generateProof() {
    try {
        generateBtn.disabled = true;
        verifyBtn.disabled = true;
        updateStatus('Generating proof...');

        const noir = new Noir(circuit as any);
        const honk = new UltraHonkBackend(noirCircuit.bytecode, { threads: 4 });

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

        console.log('Generating witness...');
        const { witness } = await noir.execute(inputs);
        
        console.log('Generating proof...');
        const result = await honk.generateProof(witness);
        
        currentProof = result.proof;
        currentPublicInputs = result.publicInputs;

        console.log('Generated proof:', {
            proofLength: currentProof.length,
            publicInputs: currentPublicInputs
        });

        updateStatus('Proof generated successfully!');
        verifyBtn.disabled = false;
    } catch (error: any) {
        console.error('Error details:', error);
        if (error instanceof Error) {
            console.error('Error stack:', error.stack);
        }
        updateStatus(`Error generating proof: ${error?.message || 'Unknown error'}`, true);
    } finally {
        generateBtn.disabled = false;
    }
}

async function verifyProof() {
    try {
        if (!currentProof || !currentPublicInputs) {
            updateStatus('Please generate a proof first', true);
            return;
        }

        verifyBtn.disabled = true;
        updateStatus('Verifying proof...');

        const honk = new UltraHonkBackend(noirCircuit.bytecode, { threads: 4 });
        
        const proofData = {
            proof: currentProof,
            publicInputs: currentPublicInputs
        };

        console.log('Verifying proof with:', {
            proofLength: proofData.proof.length,
            publicInputs: proofData.publicInputs
        });
        
        const isValid = await honk.verifyProof(proofData, { keccak: true });
        console.log('Verification result:', isValid);

        if (isValid) {
            updateStatus('Proof verified successfully!');
        } else {
            updateStatus('Proof verification failed!', true);
        }
    } catch (error: any) {
        console.error('Error details:', error);
        if (error instanceof Error) {
            console.error('Error stack:', error.stack);
        }
        updateStatus(`Error verifying proof: ${error?.message || 'Unknown error'}`, true);
    } finally {
        verifyBtn.disabled = false;
    }
}

// Make sure the status div is visible initially
statusDiv.style.display = 'block';
statusDiv.textContent = 'Ready to generate and verify proof';

generateBtn.addEventListener('click', generateProof);
verifyBtn.addEventListener('click', verifyProof);

// Initial state
verifyBtn.disabled = true;
