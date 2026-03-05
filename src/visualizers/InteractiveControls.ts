/**
 * Interactive Controls for 3D Molecular Visualization
 *
 * Phase 2: Week 7-8 - Interactive Editing
 *
 * Handles mouse interactions for atom selection, measurement,
 * and real-time editing in the Three.js molecular viewer.
 */

import * as THREE from 'three';
import { Atom, Bond } from './types';

/**
 * Selection state for atoms
 */
export interface SelectionState {
  selectedIndices: Set<number>;
  highlightedAtom: number | null;
  selectionMode: 'single' | 'multi' | 'measurement';
}

/**
 * Measurement data
 */
export interface Measurement {
  type: 'distance' | 'angle' | 'dihedral';
  atomIndices: number[];
  value: number;
  unit: string;
  label?: string;
}

/**
 * Mouse event data for raycasting
 */
export interface MouseEventData {
  x: number;
  y: number;
  button: number;
  ctrlKey: boolean;
  shiftKey: boolean;
}

/**
 * Callback types for interaction events
 */
export type AtomSelectCallback = (atomIndex: number, isSelected: boolean) => void;
export type AtomHoverCallback = (atomIndex: number | null) => void;
export type MeasurementCallback = (measurement: Measurement) => void;
export type CoordinateChangeCallback = (
  atomIndex: number,
  newPosition: { x: number; y: number; z: number }
) => void;

/**
 * Interactive controls configuration
 */
export interface InteractiveControlsConfig {
  enableSelection: boolean;
  enableMeasurement: boolean;
  enableDragToMove: boolean;
  selectionColor: string;
  highlightColor: string;
  measurementColor: string;
  hoverHighlight: boolean;
}

/**
 * Default configuration for interactive controls
 */
export const DEFAULT_INTERACTIVE_CONFIG: InteractiveControlsConfig = {
  enableSelection: true,
  enableMeasurement: true,
  enableDragToMove: false,
  selectionColor: '#FFD700', // Gold
  highlightColor: '#00FFFF', // Cyan
  measurementColor: '#FF00FF', // Magenta
  hoverHighlight: true,
};

/**
 * Interactive controls for 3D molecular visualization
 */
export class InteractiveControls {
  private raycaster: THREE.Raycaster;
  private mouse: THREE.Vector2;
  private config: InteractiveControlsConfig;
  private selectionState: SelectionState;

  // Callbacks
  private onAtomSelect?: AtomSelectCallback;
  private onAtomHover?: AtomHoverCallback;
  private onMeasurement?: MeasurementCallback;
  private onCoordinateChange?: CoordinateChangeCallback;

  // Visual indicators
  private selectionIndicators: Map<number, THREE.Mesh> = new Map();
  private measurementLines: THREE.Line[] = [];
  private measurementLabels: THREE.Sprite[] = [];

  // Drag state
  private isDragging: boolean = false;
  private draggedAtomIndex: number | null = null;
  private dragPlane: THREE.Plane | null = null;

  constructor(config?: Partial<InteractiveControlsConfig>) {
    this.raycaster = new THREE.Raycaster();
    this.mouse = new THREE.Vector2();
    this.config = { ...DEFAULT_INTERACTIVE_CONFIG, ...config };
    this.selectionState = {
      selectedIndices: new Set(),
      highlightedAtom: null,
      selectionMode: 'single',
    };
  }

  /**
   * Set callbacks for interaction events
   */
  public setCallbacks(callbacks: {
    onAtomSelect?: AtomSelectCallback;
    onAtomHover?: AtomHoverCallback;
    onMeasurement?: MeasurementCallback;
    onCoordinateChange?: CoordinateChangeCallback;
  }): void {
    this.onAtomSelect = callbacks.onAtomSelect;
    this.onAtomHover = callbacks.onAtomHover;
    this.onMeasurement = callbacks.onMeasurement;
    this.onCoordinateChange = callbacks.onCoordinateChange;
  }

  /**
   * Update mouse position from DOM event
   */
  public updateMousePosition(
    event: MouseEvent | { clientX: number; clientY: number },
    container: HTMLElement
  ): void {
    const rect = container.getBoundingClientRect();
    this.mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
    this.mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
  }

  /**
   * Get mouse position in normalized device coordinates
   */
  public getMouseNDC(): { x: number; y: number } {
    return { x: this.mouse.x, y: this.mouse.y };
  }

  /**
   * Raycast to find intersected atoms
   */
  public raycastAtoms(
    camera: THREE.Camera,
    atomMeshes: Map<number, THREE.Mesh>
  ): { atomIndex: number; distance: number } | null {
    this.raycaster.setFromCamera(this.mouse, camera);

    const meshes = Array.from(atomMeshes.values());
    const intersects = this.raycaster.intersectObjects(meshes);

    if (intersects.length > 0) {
      const intersect = intersects[0];
      // Find the atom index from the mesh
      for (const [index, mesh] of atomMeshes.entries()) {
        if (mesh === intersect.object) {
          return { atomIndex: index, distance: intersect.distance };
        }
      }
    }

    return null;
  }

  /**
   * Handle mouse click for atom selection
   */
  public handleClick(
    event: MouseEventData,
    camera: THREE.Camera,
    atomMeshes: Map<number, THREE.Mesh>,
    scene: THREE.Scene
  ): boolean {
    if (!this.config.enableSelection) {
      return false;
    }

    const result = this.raycastAtoms(camera, atomMeshes);

    if (result) {
      const { atomIndex } = result;
      const isMultiSelect = event.ctrlKey || this.selectionState.selectionMode === 'multi';

      if (this.selectionState.selectedIndices.has(atomIndex)) {
        // Deselect
        this.selectionState.selectedIndices.delete(atomIndex);
        this.removeSelectionIndicator(atomIndex, scene);
        this.onAtomSelect?.(atomIndex, false);
      } else {
        if (!isMultiSelect) {
          // Clear previous selection in single mode
          this.clearSelection(scene);
        }
        // Select
        this.selectionState.selectedIndices.add(atomIndex);
        this.addSelectionIndicator(atomIndex, atomMeshes, scene);
        this.onAtomSelect?.(atomIndex, true);
      }

      // Check if we should perform a measurement
      if (this.config.enableMeasurement) {
        this.checkAndPerformMeasurement(scene);
      }

      return true;
    } else if (!event.ctrlKey) {
      // Clicked on empty space - clear selection
      this.clearSelection(scene);
      return true;
    }

    return false;
  }

  /**
   * Handle mouse move for hover effects
   */
  public handleMouseMove(
    event: MouseEventData,
    camera: THREE.Camera,
    atomMeshes: Map<number, THREE.Mesh>,
    scene: THREE.Scene
  ): boolean {
    if (!this.config.hoverHighlight) {
      return false;
    }

    // Handle drag if active
    if (this.isDragging && this.draggedAtomIndex !== null && this.onCoordinateChange) {
      return this.handleDrag(camera, scene);
    }

    const result = this.raycastAtoms(camera, atomMeshes);

    if (result) {
      const { atomIndex } = result;
      if (this.selectionState.highlightedAtom !== atomIndex) {
        // Remove previous highlight
        if (this.selectionState.highlightedAtom !== null) {
          this.removeHighlightIndicator(this.selectionState.highlightedAtom, scene);
        }

        // Add new highlight (if not already selected)
        if (!this.selectionState.selectedIndices.has(atomIndex)) {
          this.addHighlightIndicator(atomIndex, atomMeshes, scene);
        }

        this.selectionState.highlightedAtom = atomIndex;
        this.onAtomHover?.(atomIndex);
      }
      return true;
    } else {
      // No intersection - clear highlight
      if (this.selectionState.highlightedAtom !== null) {
        this.removeHighlightIndicator(this.selectionState.highlightedAtom, scene);
        this.selectionState.highlightedAtom = null;
        this.onAtomHover?.(null);
      }
      return false;
    }
  }

  /**
   * Handle mouse down for drag operations
   */
  public handleMouseDown(
    event: MouseEventData,
    camera: THREE.Camera,
    atomMeshes: Map<number, THREE.Mesh>
  ): boolean {
    if (!this.config.enableDragToMove) {
      return false;
    }

    const result = this.raycastAtoms(camera, atomMeshes);
    if (result && this.selectionState.selectedIndices.has(result.atomIndex)) {
      this.isDragging = true;
      this.draggedAtomIndex = result.atomIndex;

      // Create drag plane facing camera
      const mesh = atomMeshes.get(result.atomIndex);
      if (mesh) {
        const normal = new THREE.Vector3();
        camera.getWorldDirection(normal);
        this.dragPlane = new THREE.Plane().setFromNormalAndCoplanarPoint(normal, mesh.position);
      }

      return true;
    }

    return false;
  }

  /**
   * Handle mouse up to end drag
   */
  public handleMouseUp(): void {
    this.isDragging = false;
    this.draggedAtomIndex = null;
    this.dragPlane = null;
  }

  /**
   * Handle drag operation
   */
  private handleDrag(camera: THREE.Camera, scene: THREE.Scene): boolean {
    if (!this.dragPlane || this.draggedAtomIndex === null) {
      return false;
    }

    this.raycaster.setFromCamera(this.mouse, camera);
    const target = new THREE.Vector3();
    this.raycaster.ray.intersectPlane(this.dragPlane, target);

    if (target) {
      this.onCoordinateChange?.(this.draggedAtomIndex, {
        x: target.x,
        y: target.y,
        z: target.z,
      });
      return true;
    }

    return false;
  }

  /**
   * Add visual indicator for selected atom
   */
  private addSelectionIndicator(
    atomIndex: number,
    atomMeshes: Map<number, THREE.Mesh>,
    scene: THREE.Scene
  ): void {
    const mesh = atomMeshes.get(atomIndex);
    if (!mesh) {
      return;
    }

    // Create selection ring (slightly larger sphere with wireframe)
    const geometry = new THREE.SphereGeometry(
      (mesh.geometry as THREE.SphereGeometry).parameters.radius * 1.2,
      32,
      32
    );
    const material = new THREE.MeshBasicMaterial({
      color: this.config.selectionColor,
      wireframe: true,
      transparent: true,
      opacity: 0.5,
    });
    const indicator = new THREE.Mesh(geometry, material);
    indicator.position.copy(mesh.position);

    scene.add(indicator);
    this.selectionIndicators.set(atomIndex, indicator);
  }

  /**
   * Remove selection indicator
   */
  private removeSelectionIndicator(atomIndex: number, scene: THREE.Scene): void {
    const indicator = this.selectionIndicators.get(atomIndex);
    if (indicator) {
      scene.remove(indicator);
      indicator.geometry.dispose();
      (indicator.material as THREE.Material).dispose();
      this.selectionIndicators.delete(atomIndex);
    }
  }

  /**
   * Add highlight indicator for hover
   */
  private addHighlightIndicator(
    atomIndex: number,
    atomMeshes: Map<number, THREE.Mesh>,
    scene: THREE.Scene
  ): void {
    // Similar to selection but with different color and no storage needed
    // Implementation would add temporary highlight
  }

  /**
   * Remove highlight indicator
   */
  private removeHighlightIndicator(atomIndex: number, scene: THREE.Scene): void {
    // Remove temporary highlight
  }

  /**
   * Clear all selections
   */
  public clearSelection(scene: THREE.Scene): void {
    this.selectionState.selectedIndices.forEach(index => {
      this.removeSelectionIndicator(index, scene);
      this.onAtomSelect?.(index, false);
    });
    this.selectionState.selectedIndices.clear();
    this.clearMeasurements(scene);
  }

  /**
   * Get currently selected atom indices
   */
  public getSelectedAtoms(): number[] {
    return Array.from(this.selectionState.selectedIndices);
  }

  /**
   * Set selection mode
   */
  public setSelectionMode(mode: 'single' | 'multi' | 'measurement'): void {
    this.selectionState.selectionMode = mode;
  }

  /**
   * Get selection mode
   */
  public getSelectionMode(): 'single' | 'multi' | 'measurement' {
    return this.selectionState.selectionMode;
  }

  /**
   * Check if measurements should be performed and do so
   */
  private checkAndPerformMeasurement(scene: THREE.Scene): void {
    const selected = this.getSelectedAtoms();

    if (selected.length === 2) {
      // Distance measurement
      this.measureDistance(selected[0], selected[1], scene);
    } else if (selected.length === 3) {
      // Angle measurement
      this.measureAngle(selected[0], selected[1], selected[2], scene);
    } else if (selected.length === 4) {
      // Dihedral measurement
      this.measureDihedral(selected[0], selected[1], selected[2], selected[3], scene);
    }
  }

  /**
   * Measure distance between two atoms
   */
  private measureDistance(atom1: number, atom2: number, scene: THREE.Scene): void {
    // This would be called by the renderer with actual atom positions
    // For now, emit an event that the renderer will handle
    this.onMeasurement?.({
      type: 'distance',
      atomIndices: [atom1, atom2],
      value: 0, // Will be calculated by renderer
      unit: 'A',
      label: `d(${atom1}-${atom2})`,
    });
  }

  /**
   * Measure angle between three atoms
   */
  private measureAngle(atom1: number, atom2: number, atom3: number, scene: THREE.Scene): void {
    this.onMeasurement?.({
      type: 'angle',
      atomIndices: [atom1, atom2, atom3],
      value: 0,
      unit: 'deg',
      label: `angle(${atom1}-${atom2}-${atom3})`,
    });
  }

  /**
   * Measure dihedral angle between four atoms
   */
  private measureDihedral(
    atom1: number,
    atom2: number,
    atom3: number,
    atom4: number,
    scene: THREE.Scene
  ): void {
    this.onMeasurement?.({
      type: 'dihedral',
      atomIndices: [atom1, atom2, atom3, atom4],
      value: 0,
      unit: 'deg',
      label: `dihedral(${atom1}-${atom2}-${atom3}-${atom4})`,
    });
  }

  /**
   * Clear all measurement visualizations
   */
  public clearMeasurements(scene: THREE.Scene): void {
    this.measurementLines.forEach(line => {
      scene.remove(line);
      line.geometry.dispose();
      (line.material as THREE.Material).dispose();
    });
    this.measurementLines = [];

    this.measurementLabels.forEach(label => {
      scene.remove(label);
      // Sprite materials are shared, don't dispose
    });
    this.measurementLabels = [];
  }

  /**
   * Draw measurement line between atoms
   */
  public drawMeasurementLine(
    positions: { x: number; y: number; z: number }[],
    scene: THREE.Scene,
    color: string = this.config.measurementColor
  ): void {
    if (positions.length < 2) {
      return;
    }

    const points = positions.map(p => new THREE.Vector3(p.x, p.y, p.z));
    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    const material = new THREE.LineBasicMaterial({ color, linewidth: 2 });
    const line = new THREE.Line(geometry, material);

    scene.add(line);
    this.measurementLines.push(line);
  }

  /**
   * Update configuration
   */
  public updateConfig(config: Partial<InteractiveControlsConfig>): void {
    this.config = { ...this.config, ...config };
  }

  /**
   * Get current configuration
   */
  public getConfig(): InteractiveControlsConfig {
    return { ...this.config };
  }

  /**
   * Dispose of all resources
   */
  public dispose(scene: THREE.Scene): void {
    this.clearSelection(scene);
    this.clearMeasurements(scene);

    if (this.selectionState.highlightedAtom !== null) {
      this.removeHighlightIndicator(this.selectionState.highlightedAtom, scene);
    }

    this.onAtomSelect = undefined;
    this.onAtomHover = undefined;
    this.onMeasurement = undefined;
    this.onCoordinateChange = undefined;
  }
}
