/**
 * ASE Converter Tests
 * 
 * Unit tests for ASE format conversion functionality.
 */

import * as assert from 'assert';
import * as path from 'path';
import * as fs from 'fs';
import { ASEConverter, ASEFormat } from '../../../ase/ASEConverter';

suite('ASE Converter Test Suite', () => {
    let converter: ASEConverter;
    const testFixturesPath = path.join(__dirname, '..', '..', 'fixtures');

    setup(() => {
        // Mock extension context
        const mockContext = {
            extensionPath: path.join(__dirname, '..', '..', '..', '..'),
            subscriptions: []
        } as any;
        
        converter = new ASEConverter(mockContext);
    });

    test('ASEConverter should be instantiated', () => {
        assert.ok(converter instanceof ASEConverter);
    });

    test('getSupportedFormats should return all supported formats', () => {
        const formats = converter.getSupportedFormats();
        
        assert.ok(formats.vasp);
        assert.ok(formats.cp2k);
        assert.ok(formats.qe);
        assert.ok(formats.gaussian);
        assert.ok(formats.orca);
        assert.ok(formats.nwchem);
        assert.ok(formats.gamess);
        assert.ok(formats.lammps);
        assert.ok(formats.xyz);
        assert.ok(formats.pdb);
        assert.ok(formats.cif);
        
        // Check structure
        assert.strictEqual(formats.vasp.name, 'VASP');
        assert.ok(Array.isArray(formats.vasp.extensions));
        assert.ok(formats.vasp.description.length > 0);
    });

    test('readToAtoms should handle missing files', async () => {
        const result = await converter.readToAtoms('/nonexistent/file.xyz');
        
        assert.strictEqual(result.success, false);
        assert.ok(result.error?.includes('File not found'));
    });

    test('convertFormat should handle invalid paths', async () => {
        const result = await converter.convertFormat(
            '/nonexistent/input.xyz',
            '/tmp/output.poscar'
        );
        
        assert.strictEqual(result.success, false);
        assert.ok(result.error);
    });

    test('isAvailable should check Python backend', async () => {
        const available = await converter.isAvailable();
        
        // This will be true if ASE is installed, false otherwise
        assert.strictEqual(typeof available, 'boolean');
    });
});
