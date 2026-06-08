/**
 * Camera controller for Three.js molecular visualization
 *
 * Manages camera position, rotation, zoom, and auto-centering.
 * Extracted from ThreeJsRenderer for maintainability.
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/34
 */

import * as THREE from 'three';
import { Atom, CameraState } from './types';

/**
 * Controls camera positioning for molecular visualization.
 */
export class CameraController {
  private state: CameraState;

  constructor(
    private camera: THREE.PerspectiveCamera,
    defaultState: CameraState
  ) {
    this.state = { ...defaultState };
  }

  /**
   * Get current camera state.
   */
  getState(): CameraState {
    return { ...this.state };
  }

  /**
   * Get camera zoom level.
   */
  getZoom(): number {
    return this.state.zoom;
  }

  /**
   * Update the default camera state and apply it.
   */
  reset(defaultState: CameraState): void {
    this.state = { ...defaultState };
    this.camera.position.set(this.state.position.x, this.state.position.y, this.state.position.z);
    this.camera.lookAt(this.state.target.x, this.state.target.y, this.state.target.z);
  }

  /**
   * Rotate camera around the target point.
   */
  rotate(deltaX: number, deltaY: number): void {
    const theta = (deltaX * Math.PI) / 180;
    const phi = (deltaY * Math.PI) / 180;

    const offset = new THREE.Vector3(
      this.camera.position.x - this.state.target.x,
      this.camera.position.y - this.state.target.y,
      this.camera.position.z - this.state.target.z
    );

    offset.applyAxisAngle(new THREE.Vector3(0, 1, 0), theta);

    const xAxis = new THREE.Vector3(1, 0, 0).applyAxisAngle(new THREE.Vector3(0, 1, 0), theta);
    offset.applyAxisAngle(xAxis, phi);

    this.camera.position.set(
      this.state.target.x + offset.x,
      this.state.target.y + offset.y,
      this.state.target.z + offset.z
    );

    this.camera.lookAt(this.state.target.x, this.state.target.y, this.state.target.z);
  }

  /**
   * Zoom camera by a multiplicative factor.
   */
  zoom(factor: number): void {
    this.state.zoom *= factor;

    const offset = new THREE.Vector3(
      this.camera.position.x - this.state.target.x,
      this.camera.position.y - this.state.target.y,
      this.camera.position.z - this.state.target.z
    );

    offset.multiplyScalar(factor);

    this.camera.position.set(
      this.state.target.x + offset.x,
      this.state.target.y + offset.y,
      this.state.target.z + offset.z
    );
  }

  /**
   * Pan camera by an offset.
   */
  pan(dx: number, dy: number, dz: number): void {
    this.camera.position.x += dx;
    this.camera.position.y += dy;
    this.camera.position.z += dz;

    this.state.target.x += dx;
    this.state.target.y += dy;
    this.state.target.z += dz;
  }

  /**
   * Center camera on a molecular structure defined by its atoms.
   */
  centerOnAtoms(atoms: Atom[]): void {
    if (atoms.length === 0) {
      return;
    }

    const box = new THREE.Box3();
    atoms.forEach(atom => {
      box.expandByPoint(new THREE.Vector3(atom.x, atom.y, atom.z));
    });

    const center = box.getCenter(new THREE.Vector3());

    this.state.target = { x: center.x, y: center.y, z: center.z };

    const size = box.getSize(new THREE.Vector3());
    const maxDim = Math.max(size.x, size.y, size.z);
    const distance = maxDim * 2;

    this.camera.position.set(center.x + distance, center.y + distance, center.z + distance);
    this.camera.lookAt(center);
  }
}
