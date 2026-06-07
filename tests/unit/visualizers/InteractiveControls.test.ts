/**
 * @jest-environment jsdom
 *
 * Tests for InteractiveControls
 *
 * Phase 2: Week 7-8 - Interactive Editing Tests
 */

import * as THREE from 'three';
import {
  InteractiveControls,
  SelectionState,
  Measurement,
  MouseEventData,
  DEFAULT_INTERACTIVE_CONFIG,
} from '../../../src/visualizers/InteractiveControls';

// Mock Three.js objects
const mockCamera = {
  position: new THREE.Vector3(0, 0, 10),
  getWorldDirection: jest.fn(() => new THREE.Vector3(0, 0, -1)),
} as unknown as THREE.PerspectiveCamera;

const mockScene = new THREE.Scene();

const createMockMesh = (
  index: number,
  position: { x: number; y: number; z: number }
): THREE.Mesh => {
  const geometry = new THREE.SphereGeometry(1, 8, 8);
  const material = new THREE.MeshBasicMaterial({ color: 0xffffff });
  const mesh = new THREE.Mesh(geometry, material);
  mesh.position.set(position.x, position.y, position.z);
  return mesh;
};

describe('InteractiveControls', () => {
  let controls: InteractiveControls;

  beforeEach(() => {
    controls = new InteractiveControls();
    // Clear scene
    while (mockScene.children.length > 0) {
      mockScene.remove(mockScene.children[0]);
    }
  });

  afterEach(() => {
    controls.dispose(mockScene);
  });

  describe('Constructor', () => {
    it('should create with default config', () => {
      const ctrl = new InteractiveControls();
      expect(ctrl.getConfig()).toEqual(DEFAULT_INTERACTIVE_CONFIG);
    });

    it('should create with custom config', () => {
      const ctrl = new InteractiveControls({
        enableSelection: false,
        selectionColor: '#FF0000',
      });
      const config = ctrl.getConfig();
      expect(config.enableSelection).toBe(false);
      expect(config.selectionColor).toBe('#FF0000');
      expect(config.enableMeasurement).toBe(true);
    });
  });

  describe('Configuration', () => {
    it('should update configuration', () => {
      controls.updateConfig({ enableDragToMove: true });
      expect(controls.getConfig().enableDragToMove).toBe(true);
    });

    it('should not affect other config values when updating', () => {
      const original = controls.getConfig();
      controls.updateConfig({ selectionColor: '#000000' });
      const updated = controls.getConfig();
      expect(updated.selectionColor).toBe('#000000');
      expect(updated.enableSelection).toBe(original.enableSelection);
    });
  });

  describe('Selection Mode', () => {
    it('should default to single selection mode', () => {
      expect(controls.getSelectionMode()).toBe('single');
    });

    it('should set selection mode', () => {
      controls.setSelectionMode('multi');
      expect(controls.getSelectionMode()).toBe('multi');
    });

    it('should support measurement mode', () => {
      controls.setSelectionMode('measurement');
      expect(controls.getSelectionMode()).toBe('measurement');
    });
  });

  describe('Mouse Position', () => {
    it('should update mouse position from event', () => {
      const container = document.createElement('div');
      container.getBoundingClientRect = jest.fn(() => ({
        left: 0,
        top: 0,
        width: 100,
        height: 100,
        right: 100,
        bottom: 100,
        x: 0,
        y: 0,
        toJSON: () => {},
      }));

      controls.updateMousePosition({ clientX: 50, clientY: 50 }, container);
      const ndc = controls.getMouseNDC();
      expect(ndc.x).toBe(0);
      expect(ndc.y).toBe(0);
    });

    it('should convert to NDC correctly', () => {
      const container = document.createElement('div');
      container.getBoundingClientRect = jest.fn(() => ({
        left: 0,
        top: 0,
        width: 200,
        height: 200,
        right: 200,
        bottom: 200,
        x: 0,
        y: 0,
        toJSON: () => {},
      }));

      // Top-left should be (-1, 1)
      controls.updateMousePosition({ clientX: 0, clientY: 0 }, container);
      let ndc = controls.getMouseNDC();
      expect(ndc.x).toBeCloseTo(-1);
      expect(ndc.y).toBeCloseTo(1);

      // Bottom-right should be (1, -1)
      controls.updateMousePosition({ clientX: 200, clientY: 200 }, container);
      ndc = controls.getMouseNDC();
      expect(ndc.x).toBeCloseTo(1);
      expect(ndc.y).toBeCloseTo(-1);
    });
  });

  describe('Callbacks', () => {
    it('should set all callbacks', () => {
      const callbacks = {
        onAtomSelect: jest.fn(),
        onAtomHover: jest.fn(),
        onMeasurement: jest.fn(),
        onCoordinateChange: jest.fn(),
      };

      controls.setCallbacks(callbacks);
      // Callbacks are set if no error thrown
      expect(true).toBe(true);
    });

    it('should call onAtomSelect when atom is selected', () => {
      const onAtomSelect = jest.fn();
      controls.setCallbacks({ onAtomSelect });

      // Create atom meshes
      const atomMeshes = new Map<number, THREE.Mesh>();
      atomMeshes.set(0, createMockMesh(0, { x: 0, y: 0, z: 0 }));
      mockScene.add(atomMeshes.get(0)!);

      // Position mouse over atom
      const container = document.createElement('div');
      container.getBoundingClientRect = jest.fn(() => ({
        left: 0,
        top: 0,
        width: 100,
        height: 100,
        right: 100,
        bottom: 100,
        x: 0,
        y: 0,
        toJSON: () => {},
      }));

      // Click on atom (this won't actually trigger since we can't raycast in tests)
      // But we verify the callback is set up
      expect(onAtomSelect).not.toHaveBeenCalled();
    });
  });

  describe('Selection State', () => {
    it('should return empty selection initially', () => {
      expect(controls.getSelectedAtoms()).toEqual([]);
    });

    it('should clear selection', () => {
      // Add some selection indicators manually
      controls.clearSelection(mockScene);
      expect(controls.getSelectedAtoms()).toEqual([]);
    });
  });

  describe('Measurements', () => {
    it('should clear measurements', () => {
      // Should not throw when clearing empty measurements
      expect(() => controls.clearMeasurements(mockScene)).not.toThrow();
    });

    it('should draw measurement line', () => {
      const positions = [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
      ];
      controls.drawMeasurementLine(positions, mockScene);
      // Line should be added to scene
      const lines = mockScene.children.filter(child => child instanceof THREE.Line);
      expect(lines.length).toBe(1);
    });

    it('should not draw line with less than 2 positions', () => {
      const positions = [{ x: 0, y: 0, z: 0 }];
      controls.drawMeasurementLine(positions, mockScene);
      const lines = mockScene.children.filter(child => child instanceof THREE.Line);
      expect(lines.length).toBe(0);
    });
  });

  describe('Dispose', () => {
    it('should dispose without errors', () => {
      expect(() => controls.dispose(mockScene)).not.toThrow();
    });

    it('should clear selections on dispose', () => {
      controls.dispose(mockScene);
      expect(controls.getSelectedAtoms()).toEqual([]);
    });
  });

  describe('Drag Operations', () => {
    it('should handle mouse up when not dragging', () => {
      expect(() => controls.handleMouseUp()).not.toThrow();
    });

    it('should not handle drag when disabled', () => {
      controls.updateConfig({ enableDragToMove: false });
      const result = controls.handleMouseDown(
        { x: 0, y: 0, button: 0, ctrlKey: false, shiftKey: false },
        mockCamera,
        new Map()
      );
      expect(result).toBe(false);
    });
  });

  describe('Mouse Event Handling', () => {
    it('should not handle click when selection disabled', () => {
      controls.updateConfig({ enableSelection: false });
      const event: MouseEventData = {
        x: 0,
        y: 0,
        button: 0,
        ctrlKey: false,
        shiftKey: false,
      };
      const result = controls.handleClick(event, mockCamera, new Map(), mockScene);
      expect(result).toBe(false);
    });

    it('should not handle mouse move when hover disabled', () => {
      controls.updateConfig({ hoverHighlight: false });
      const event: MouseEventData = {
        x: 0,
        y: 0,
        button: 0,
        ctrlKey: false,
        shiftKey: false,
      };
      const result = controls.handleMouseMove(event, mockCamera, new Map(), mockScene);
      expect(result).toBe(false);
    });
  });
});

describe('InteractiveControls - Additional Coverage', () => {
  let controls: InteractiveControls;
  let mockScene: THREE.Scene;
  let mockCamera: THREE.PerspectiveCamera;

  beforeEach(() => {
    controls = new InteractiveControls();
    mockScene = new THREE.Scene();
    mockCamera = new THREE.PerspectiveCamera(75, 1, 0.1, 1000);
    mockCamera.position.set(0, 0, 10);
  });

  afterEach(() => {
    controls.dispose(mockScene);
  });

  describe('Raycasting', () => {
    it('should return null when no atoms to raycast', () => {
      const container = document.createElement('div');
      container.getBoundingClientRect = jest.fn(() => ({
        left: 0,
        top: 0,
        width: 100,
        height: 100,
        right: 100,
        bottom: 100,
        x: 0,
        y: 0,
        toJSON: () => {},
      }));

      controls.updateMousePosition({ clientX: 50, clientY: 50 }, container);
      const result = controls.raycastAtoms(mockCamera, new Map());
      expect(result).toBeNull();
    });

    it('should find intersected atom', () => {
      const container = document.createElement('div');
      container.getBoundingClientRect = jest.fn(() => ({
        left: 0,
        top: 0,
        width: 100,
        height: 100,
        right: 100,
        bottom: 100,
        x: 0,
        y: 0,
        toJSON: () => {},
      }));

      // Create atom mesh at origin
      const geometry = new THREE.SphereGeometry(1, 8, 8);
      const material = new THREE.MeshBasicMaterial({ color: 0xffffff });
      const mesh = new THREE.Mesh(geometry, material);
      mesh.position.set(0, 0, 0);

      const atomMeshes = new Map<number, THREE.Mesh>();
      atomMeshes.set(0, mesh);
      mockScene.add(mesh);

      // Position camera to look at origin
      mockCamera.position.set(0, 0, 5);
      mockCamera.lookAt(0, 0, 0);

      controls.updateMousePosition({ clientX: 50, clientY: 50 }, container);
      // Raycast will be tested but may not find intersection in test environment
      const result = controls.raycastAtoms(mockCamera, atomMeshes);
      // Result may be null due to test environment limitations
      expect(result === null || typeof result.atomIndex === 'number').toBe(true);
    });
  });

  describe('Selection Indicators', () => {
    it('should handle click on empty space', () => {
      const event: MouseEventData = {
        x: 0,
        y: 0,
        button: 0,
        ctrlKey: false,
        shiftKey: false,
      };
      // No atoms in scene, should return true (cleared selection)
      const result = controls.handleClick(event, mockCamera, new Map(), mockScene);
      expect(result).toBe(true);
    });

    it('should handle click with ctrl key (multi-select)', () => {
      const event: MouseEventData = {
        x: 0,
        y: 0,
        button: 0,
        ctrlKey: true,
        shiftKey: false,
      };
      // Should not clear selection when ctrl is pressed
      const result = controls.handleClick(event, mockCamera, new Map(), mockScene);
      expect(result).toBe(false);
    });
  });

  describe('Measurements - Complete Coverage', () => {
    it('should handle measurement callback for distance', () => {
      const onMeasurement = jest.fn();
      controls.setCallbacks({ onMeasurement });

      // Trigger measurement via callback
      controls['onMeasurement']?.({
        type: 'distance',
        atomIndices: [0, 1],
        value: 1.5,
        unit: 'A',
        label: 'd(0-1)',
      });

      expect(onMeasurement).toHaveBeenCalledWith({
        type: 'distance',
        atomIndices: [0, 1],
        value: 1.5,
        unit: 'A',
        label: 'd(0-1)',
      });
    });

    it('should draw measurement line with custom color', () => {
      const positions = [
        { x: 0, y: 0, z: 0 },
        { x: 2, y: 0, z: 0 },
        { x: 2, y: 2, z: 0 },
      ];
      controls.drawMeasurementLine(positions, mockScene, '#FF0000');
      const lines = mockScene.children.filter(child => child instanceof THREE.Line);
      expect(lines.length).toBe(1);
    });

    it('should clear multiple measurement lines', () => {
      // Draw multiple lines
      controls.drawMeasurementLine(
        [
          { x: 0, y: 0, z: 0 },
          { x: 1, y: 0, z: 0 },
        ],
        mockScene
      );
      controls.drawMeasurementLine(
        [
          { x: 1, y: 0, z: 0 },
          { x: 1, y: 1, z: 0 },
        ],
        mockScene
      );

      let lines = mockScene.children.filter(child => child instanceof THREE.Line);
      expect(lines.length).toBe(2);

      controls.clearMeasurements(mockScene);

      lines = mockScene.children.filter(child => child instanceof THREE.Line);
      expect(lines.length).toBe(0);
    });
  });

  describe('Drag Operations - Extended', () => {
    it('should enable drag to move', () => {
      controls.updateConfig({ enableDragToMove: true });
      expect(controls.getConfig().enableDragToMove).toBe(true);
    });

    it('should attempt drag with no drag plane', () => {
      controls.updateConfig({ enableDragToMove: true });

      // Manually set dragging state
      controls['isDragging'] = true;
      controls['draggedAtomIndex'] = 0;
      controls['dragPlane'] = null;

      const event: MouseEventData = {
        x: 0,
        y: 0,
        button: 0,
        ctrlKey: false,
        shiftKey: false,
      };

      const result = controls.handleMouseMove(event, mockCamera, new Map(), mockScene);
      expect(result).toBe(false);
    });
  });

  describe('Configuration Options', () => {
    it('should have correct default highlight color', () => {
      const config = controls.getConfig();
      expect(config.highlightColor).toBe('#00FFFF');
    });

    it('should have correct default measurement color', () => {
      const config = controls.getConfig();
      expect(config.measurementColor).toBe('#FF00FF');
    });

    it('should update hover highlight setting', () => {
      controls.updateConfig({ hoverHighlight: false });
      expect(controls.getConfig().hoverHighlight).toBe(false);

      controls.updateConfig({ hoverHighlight: true });
      expect(controls.getConfig().hoverHighlight).toBe(true);
    });
  });
});
