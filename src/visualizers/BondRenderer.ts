/**
 * Bond rendering for Three.js molecular visualization
 *
 * Handles bond detection, mesh creation, and disposal.
 * Extracted from ThreeJsRenderer for maintainability.
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/34
 */

import * as THREE from 'three';
import { Atom, Bond, COVALENT_RADII } from './types';

/**
 * Bond detection tolerance factor.
 * Bonds are detected when interatomic distance is within
 * this factor of the sum of covalent radii.
 */
const BOND_TOLERANCE = 1.2;

/**
 * Manages bond mesh lifecycle within a Three.js scene.
 */
export class BondRenderer {
  private bondMeshes: Map<string, THREE.Mesh> = new Map();

  constructor(
    private scene: THREE.Scene,
    private bondRadius: number
  ) {}

  /**
   * Update the bond radius config.
   */
  updateBondRadius(radius: number): void {
    this.bondRadius = radius;
  }

  /**
   * Detect and render bonds between atoms.
   */
  detectAndRenderBonds(atoms: Atom[], showBonds: boolean): Bond[] {
    if (!showBonds) {
      return [];
    }

    this.clearBonds();

    const bonds = this.detectBonds(atoms);

    bonds.forEach(bond => {
      const mesh = this.createBondMesh(bond, atoms);
      if (mesh) {
        this.bondMeshes.set(`${bond.atomIndex1}-${bond.atomIndex2}`, mesh);
        this.scene.add(mesh);
      }
    });

    return bonds;
  }

  /**
   * Get the number of currently rendered bonds.
   */
  getBondCount(): number {
    return this.bondMeshes.size;
  }

  /**
   * Get bond threshold distance for two elements.
   */
  static getBondThreshold(element1: string, element2: string): number {
    return (
      (COVALENT_RADII[element1] || COVALENT_RADII.default) +
      (COVALENT_RADII[element2] || COVALENT_RADII.default)
    );
  }

  /**
   * Clear all bond meshes from the scene and dispose resources.
   */
  clearBonds(): void {
    this.bondMeshes.forEach(mesh => {
      this.scene.remove(mesh);
      mesh.geometry.dispose();
      if (Array.isArray(mesh.material)) {
        mesh.material.forEach(m => m.dispose());
      } else {
        mesh.material.dispose();
      }
    });
    this.bondMeshes.clear();
  }

  /**
   * Detect bonds between atoms based on distance and covalent radii.
   */
  private detectBonds(atoms: Atom[]): Bond[] {
    const bonds: Bond[] = [];

    for (let i = 0; i < atoms.length; i++) {
      for (let j = i + 1; j < atoms.length; j++) {
        const atom1 = atoms[i];
        const atom2 = atoms[j];

        const distance = this.calculateDistance(atom1, atom2);
        const maxDistance = BondRenderer.getBondThreshold(atom1.element, atom2.element);

        if (distance <= maxDistance * BOND_TOLERANCE) {
          bonds.push({
            atomIndex1: i,
            atomIndex2: j,
            length: distance,
          });
        }
      }
    }

    return bonds;
  }

  /**
   * Calculate Euclidean distance between two atoms.
   */
  private calculateDistance(atom1: Atom, atom2: Atom): number {
    const dx = atom2.x - atom1.x;
    const dy = atom2.y - atom1.y;
    const dz = atom2.z - atom1.z;

    return Math.sqrt(dx * dx + dy * dy + dz * dz);
  }

  /**
   * Create a cylinder mesh for a bond.
   */
  private createBondMesh(bond: Bond, atoms: Atom[]): THREE.Mesh | null {
    const atom1 = atoms[bond.atomIndex1];
    const atom2 = atoms[bond.atomIndex2];

    if (!atom1 || !atom2) {
      return null;
    }

    const direction = new THREE.Vector3(atom2.x - atom1.x, atom2.y - atom1.y, atom2.z - atom1.z);
    const length = direction.length();

    const geometry = new THREE.CylinderGeometry(this.bondRadius, this.bondRadius, length, 8);

    const material = new THREE.MeshPhongMaterial({
      color: 0x888888,
    });

    const mesh = new THREE.Mesh(geometry, material);

    const midpoint = new THREE.Vector3(
      (atom1.x + atom2.x) / 2,
      (atom1.y + atom2.y) / 2,
      (atom1.z + atom2.z) / 2
    );
    mesh.position.copy(midpoint);
    mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.normalize());

    return mesh;
  }
}
