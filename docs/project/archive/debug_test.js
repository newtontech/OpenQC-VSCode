const path = require('path');
const fs = require('fs');
const { ASEConverter, ASEFormat } = require('./out/ase/ASEConverter');

// Mock extension context
class MockExtensionContext {
  constructor(extensionPath) {
    this.extensionPath = extensionPath;
    this.subscriptions = [];
  }
  asAbsolutePath(relativePath) {
    return path.join(this.extensionPath, relativePath);
  }
}

async function test() {
  const mockContext = new MockExtensionContext(__dirname);
  const converter = new ASEConverter(mockContext);
  
  // Check if backend is available
  const available = await converter.isAvailable();
  console.log('Backend available:', available);
  
  if (!available) {
    console.log('Backend not available, exiting');
    return;
  }
  
  // Test conversion
  const inputFile = path.join(__dirname, 'tests/fixtures/migration/POSCAR');
  const outputFile = path.join(__dirname, 'tests/temp/migration/debug_test.xyz');
  
  // Ensure output directory exists
  const outputDir = path.dirname(outputFile);
  if (!fs.existsSync(outputDir)) {
    fs.mkdirSync(outputDir, { recursive: true });
  }
  
  console.log('Input file:', inputFile);
  console.log('Output file:', outputFile);
  console.log('Input exists:', fs.existsSync(inputFile));
  
  const result = await converter.convertFormat(
    inputFile,
    outputFile,
    ASEFormat.VASP,
    ASEFormat.XYZ
  );
  
  console.log('Conversion result:', JSON.stringify(result, null, 2));
}

test().catch(console.error);
