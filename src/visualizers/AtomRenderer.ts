/**
 * Atom rendering for Three.js molecular visualization
 *
 * Handles creation, storage, and disposal of atom meshes.
 * Extracted from ThreeJsRenderer for maintainability.
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/34
 */

import * as THREE from 'three';
import { Atom, VisualizationConfig, ELEMENT_COLORS, COVALENT_RADII, VDW_RADII } from './types';

/**
 * Manages atom mesh lifecycle within a Three.js scene.
 */
export class AtomRenderer {
  private atomMeshes: Map<number, THREE.Mesh> = new Map();

  constructor(
    private scene: THREE.Scene,
    private config: VisualizationConfig
  ) {}

  /**
   * Update the visualization config reference (e.g. after mode change).
   */
  updateConfig(config: VisualizationConfig): void {
    this.config = config;
  }

  /**
   * Render all atoms to the scene, returning render timing info.
   */
  renderAtoms(atoms: Atom[]): { atomCount: number; renderTime: number } {
    const startTime = performance.now();

    this.clearAtoms();

    atoms.forEach((atom, index) => {
      const mesh = this.createAtomMesh(atom);
      if (mesh) {
        this.atomMeshes.set(index, mesh);
        this.scene.add(mesh);
      }
    });

    const renderTime = performance.now() - startTime;
    return { atomCount: atoms.length, renderTime };
  }

  /**
   * Get the number of currently rendered atoms.
   */
  getAtomCount(): number {
    return this.atomMeshes.size;
  }

  /**
   * Get current atom colors for all unique elements.
   */
  getAtomColors(atoms: Atom[]): Array<{ element: string; color: string }> {
    const uniqueElements = [...new Set(atoms.map(a => a.element))];
    return uniqueElements.map(element => ({
      element,
      color: this.getAtomColor(element),
    }));
  }

  /**
   * Get current atom radii for all unique elements.
   */
  getAtomRadii(atoms: Atom[]): Record<string, number> {
    const uniqueElements = [...new Set(atoms.map(a => a.element))];
    const radii: Record<string, number> = {};
    uniqueElements.forEach(element => {
      radii[element] = this.getAtomRadius(element);
    });
    return radii;
  }

  /**
   * Clear all atom meshes from the scene and dispose resources.
   */
  clearAtoms(): void {
    this.atomMeshes.forEach(mesh => {
      this.scene.remove(mesh);
      mesh.geometry.dispose();
      if (Array.isArray(mesh.material)) {
        mesh.material.forEach(m => m.dispose());
      } else {
        mesh.material.dispose();
      }
    });
    this.atomMeshes.clear();
  }

  /**
   * Create a mesh for a single atom.
   */
  private createAtomMesh(atom: Atom): THREE.Mesh | null {
    const radius = this.getAtomRadius(atom.element);
    const color = this.getAtomColor(atom.element);

    const geometry = new THREE.SphereGeometry(radius, 32, 32);
    const material = new THREE.MeshPhongMaterial({
      color,
      shininess: 100,
    });

    const mesh = new THREE.Mesh(geometry, material);
    mesh.position.set(atom.x, atom.y, atom.z);

    return mesh;
  }

  /**
   * Get radius for an atom based on representation mode.
   */
  private getAtomRadius(element: string): number {
    const baseRadius =
      this.config.representationMode === 'space-filling'
        ? VDW_RADII[element] || VDW_RADII.default
        : COVALENT_RADII[element] || COVALENT_RADII.default;

    return baseRadius * this.config.atomScale;
  }

  /**
   * Get color for an element.
   */
  private getAtomColor(element: string): string {
    return ELEMENT_COLORS[element] || ELEMENT_COLORS.default;
  }
}
