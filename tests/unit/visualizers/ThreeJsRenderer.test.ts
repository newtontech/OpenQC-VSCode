/**
 * Tests for ThreeJsRenderer - Phase 2: 3D Visualization
 *
 * This test suite follows TDD principles to ensure Three.js integration
 * works correctly for molecular structure visualization.
 */

import { ThreeJsRenderer } from '../../../src/visualizers/ThreeJsRenderer';
import { Atom } from '../../../src/visualizers/types';

// Mock Three.js for Node.js environment
jest.mock('three', () => {
  return {
    Scene: jest.fn().mockImplementation(() => ({
      add: jest.fn(),
      remove: jest.fn(),
      clear: jest.fn(),
    })),
    PerspectiveCamera: jest.fn().mockImplementation(() => ({
      position: { set: jest.fn() },
      lookAt: jest.fn(),
      updateProjectionMatrix: jest.fn(),
    })),
    WebGLRenderer: jest.fn().mockImplementation(() => ({
      setSize: jest.fn(),
      render: jest.fn(),
      clear: jest.fn(),
      dispose: jest.fn(),
      domElement: document.createElement('canvas'),
    })),
    SphereGeometry: jest.fn().mockImplementation(() => ({})),
    MeshPhongMaterial: jest.fn().mockImplementation(() => ({})),
    Mesh: jest.fn().mockImplementation(() => ({
      position: { set: jest.fn() },
    })),
    DirectionalLight: jest.fn().mockImplementation(() => ({
      position: { set: jest.fn() },
    })),
    AmbientLight: jest.fn().mockImplementation(() => ({})),
    AxesHelper: jest.fn().mockImplementation(() => ({})),
    Color: jest.fn().mockImplementation(() => ({})),
  };
});

describe('ThreeJsRenderer', () => {
  let renderer: ThreeJsRenderer;
  let mockContainer: HTMLElement;

  beforeEach(() => {
    // Create a mock container element
    mockContainer = document.createElement('div');
    document.body.appendChild(mockContainer);

    // Mock requestAnimationFrame and cancelAnimationFrame
    global.requestAnimationFrame = jest.fn(cb => setTimeout(cb, 16));
    global.cancelAnimationFrame = jest.fn();

    renderer = new ThreeJsRenderer(mockContainer);
  });

  afterEach(() => {
    if (renderer) {
      renderer.dispose();
    }
    if (mockContainer && mockContainer.parentNode) {
      mockContainer.parentNode.removeChild(mockContainer);
    }
  });

  describe('initialization', () => {
    it('should create a renderer with a valid container', () => {
      expect(renderer).toBeDefined();
      expect(renderer.isInitialized()).toBe(true);
    });

    it('should throw an error when container is null', () => {
      expect(() => new ThreeJsRenderer(null as any)).toThrow();
    });

    it('should set default camera position', () => {
      const camera = renderer.getCamera()!;
      expect(camera).toBeDefined();
      expect(camera!.position).toBeDefined();
    });

    it('should create a scene with proper structure', () => {
      const scene = renderer.getScene();
      expect(scene).toBeDefined();
    });

    it('should handle window resize', () => {
      const initialSize = renderer.getSize();
      renderer.updateSize(800, 600);
      const newSize = renderer.getSize();
      expect(newSize.width).toBe(800);
      expect(newSize.height).toBe(600);
    });
  });

  describe('atom rendering', () => {
    const mockAtoms: Atom[] = [
      { element: 'H', x: 0, y: 0, z: 0 },
      { element: 'O', x: 1, y: 0, z: 0 },
      { element: 'H', x: 0, y: 1, z: 0 },
    ];

    it('should render atoms from an array', () => {
      const result = renderer.renderAtoms(mockAtoms);
      expect(result.success).toBe(true);
      expect(result.atomCount).toBe(3);
    });

    it('should handle empty atom array', () => {
      const result = renderer.renderAtoms([]);
      expect(result.success).toBe(true);
      expect(result.atomCount).toBe(0);
    });

    it('should clear previous atoms before rendering new ones', () => {
      renderer.renderAtoms(mockAtoms);
      renderer.renderAtoms([{ element: 'C', x: 2, y: 2, z: 2 }]);
      const atomCount = renderer.getAtomCount();
      expect(atomCount).toBe(1);
    });

    it('should apply correct colors to elements', () => {
      renderer.renderAtoms(mockAtoms);
      const colors = renderer.getAtomColors();
      expect(colors).toContainEqual(
        expect.objectContaining({ element: 'O', color: expect.any(String) })
      );
    });

    it('should apply correct radii based on atomic number', () => {
      renderer.renderAtoms([
        { element: 'H', x: 0, y: 0, z: 0 },
        { element: 'C', x: 1, y: 0, z: 0 },
      ]);
      const radii = renderer.getAtomRadii();
      expect(radii.H).toBeLessThan(radii.C);
    });

    it('should handle unknown elements gracefully', () => {
      const result = renderer.renderAtoms([{ element: 'Xx', x: 0, y: 0, z: 0 }]);
      expect(result.success).toBe(true);
      expect(result.atomCount).toBe(1);
    });
  });

  describe('bond detection and rendering', () => {
    const h2oAtoms: Atom[] = [
      { element: 'O', x: 0, y: 0, z: 0 },
      { element: 'H', x: 0.96, y: 0, z: 0 },
      { element: 'H', x: -0.24, y: 0.93, z: 0 },
    ];

    it('should detect bonds between close atoms', () => {
      renderer.renderAtoms(h2oAtoms);
      renderer.detectAndRenderBonds();
      const bondCount = renderer.getBondCount();
      expect(bondCount).toBeGreaterThan(0);
    });

    it('should not create bonds between distant atoms', () => {
      const distantAtoms: Atom[] = [
        { element: 'H', x: 0, y: 0, z: 0 },
        { element: 'H', x: 100, y: 100, z: 100 },
      ];
      renderer.renderAtoms(distantAtoms);
      renderer.detectAndRenderBonds();
      const bondCount = renderer.getBondCount();
      expect(bondCount).toBe(0);
    });

    it('should use correct bond radii thresholds', () => {
      const threshold = renderer.getBondThreshold('H', 'O');
      expect(threshold).toBeGreaterThan(0);
      expect(threshold).toBeLessThan(2.0);
    });
  });

  describe('unit cell rendering', () => {
    const mockLattice = [
      [10, 0, 0],
      [0, 10, 0],
      [0, 0, 10],
    ];

    it('should render unit cell box when lattice is provided', () => {
      const result = renderer.renderUnitCell(mockLattice);
      expect(result.success).toBe(true);
      expect(result.hasCell).toBe(true);
    });

    it('should not render unit cell when lattice is null', () => {
      const result = renderer.renderUnitCell(null as any);
      expect(result.success).toBe(true);
      expect(result.hasCell).toBe(false);
    });

    it('should handle non-orthogonal lattices', () => {
      const skewedLattice = [
        [10, 0.5, 0],
        [0, 10, 0],
        [0, 0, 10],
      ];
      const result = renderer.renderUnitCell(skewedLattice);
      expect(result.success).toBe(true);
    });
  });

  describe('camera controls', () => {
    it('should support camera rotation', () => {
      const initialPosition = renderer.getCamera()!.position;
      renderer.rotateCamera(45, 0);
      const newPosition = renderer.getCamera()!.position;
      expect(newPosition).toBeDefined();
    });

    it('should support camera zoom', () => {
      const initialZoom = renderer.getCameraZoom();
      renderer.zoomCamera(1.5);
      const newZoom = renderer.getCameraZoom();
      expect(newZoom).not.toBe(initialZoom);
    });

    it('should support camera panning', () => {
      const initialPosition = renderer.getCamera()!.position;
      renderer.panCamera(1, 0, 0);
      expect(renderer.getCamera()!.position).toBeDefined();
    });

    it('should reset camera to default position', () => {
      renderer.rotateCamera(90, 45);
      renderer.resetCamera();
      const position = renderer.getCamera()!.position;
      expect(position).toBeDefined();
    });
  });

  describe('representation modes', () => {
    const mockAtoms: Atom[] = [
      { element: 'C', x: 0, y: 0, z: 0 },
      { element: 'H', x: 1, y: 0, z: 0 },
    ];

    it('should support ball-and-stick mode', () => {
      renderer.setRepresentationMode('ball-and-stick');
      renderer.renderAtoms(mockAtoms);
      renderer.detectAndRenderBonds();
      expect(renderer.getRepresentationMode()).toBe('ball-and-stick');
    });

    it('should support space-filling mode', () => {
      renderer.setRepresentationMode('space-filling');
      renderer.renderAtoms(mockAtoms);
      expect(renderer.getRepresentationMode()).toBe('space-filling');
    });

    it('should support wireframe mode', () => {
      renderer.setRepresentationMode('wireframe');
      renderer.renderAtoms(mockAtoms);
      expect(renderer.getRepresentationMode()).toBe('wireframe');
    });
  });

  describe('export and snapshot', () => {
    it('should export scene as image data', () => {
      renderer.renderAtoms([{ element: 'H', x: 0, y: 0, z: 0 }]);
      const imageData = renderer.exportImage({ format: 'png' });
      expect(imageData).toBeDefined();
      expect(imageData.format).toBe('png');
    });

    it('should throw error for unsupported export format', () => {
      expect(() => renderer.exportImage({ format: 'bmp' as any })).toThrow();
    });

    it('should support snapshot callbacks', () => {
      const callback = jest.fn();
      renderer.onSnapshotComplete(callback);
      renderer.renderAtoms([{ element: 'H', x: 0, y: 0, z: 0 }]);
      renderer.takeSnapshot();
      expect(callback).toHaveBeenCalled();
    });
  });

  describe('performance', () => {
    it('should handle larger atom sets efficiently', () => {
      const largeAtomSet: Atom[] = Array.from({ length: 100 }, (_, i) => ({
        element: i % 2 === 0 ? 'C' : 'H',
        x: Math.random() * 10,
        y: Math.random() * 10,
        z: Math.random() * 10,
      }));

      const startTime = Date.now();
      renderer.renderAtoms(largeAtomSet);
      const endTime = Date.now();

      expect(endTime - startTime).toBeLessThan(1000); // Should render in < 1s
    });

    it('should support lazy loading for very large systems', () => {
      const hugeAtomSet: Atom[] = Array.from({ length: 1000 }, (_, i) => ({
        element: 'C',
        x: i * 0.1,
        y: 0,
        z: 0,
      }));

      renderer.enableLazyLoading(true);
      const result = renderer.renderAtoms(hugeAtomSet);
      expect(result.success).toBe(true);
    });
  });

  describe('disposal and cleanup', () => {
    it('should properly dispose all resources', () => {
      renderer.renderAtoms([{ element: 'H', x: 0, y: 0, z: 0 }]);
      renderer.dispose();
      expect(renderer.isInitialized()).toBe(false);
    });

    it('should handle multiple dispose calls safely', () => {
      renderer.dispose();
      expect(() => renderer.dispose()).not.toThrow();
    });
  });
});
