/**
 * Three.js-based 3D molecular structure renderer
 *
 * Phase 2: 3D Visualization - Core Renderer
 *
 * This module provides Three.js integration for rendering molecular structures
 * from quantum chemistry input files (VASP POSCAR, Gaussian, ORCA, etc.)
 *
 * Internal concerns are delegated to focused sub-modules:
 * - AtomRenderer: atom mesh creation and lifecycle
 * - BondRenderer: bond detection and mesh lifecycle
 * - CameraController: camera positioning and movement
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/34
 */

// DOM type declarations for Node.js environment
/// <reference lib="dom" />

import * as THREE from 'three';
import {
  Atom,
  Bond,
  MolecularStructure,
  RepresentationMode,
  CameraState,
  RenderResult,
  ExportOptions,
  VisualizationConfig,
  RendererState,
} from './types';
import { AtomRenderer } from './AtomRenderer';
import { BondRenderer } from './BondRenderer';
import { CameraController } from './CameraController';

/**
 * Three.js-based molecular structure renderer
 */
export class ThreeJsRenderer {
  private scene: THREE.Scene | null = null;
  private camera: THREE.PerspectiveCamera | null = null;
  private renderer: THREE.WebGLRenderer | null = null;
  private container: HTMLElement | null = null;
  private animationId: number | null = null;

  // Sub-renderers
  private atomRenderer: AtomRenderer | null = null;
  private bondRenderer: BondRenderer | null = null;
  private cameraController: CameraController | null = null;

  // 3D objects
  private unitCellLines: THREE.LineSegments | null = null;
  private axesHelper: THREE.AxesHelper | null = null;

  // State
  private atoms: Atom[] = [];
  private bonds: Bond[] = [];
  private config: VisualizationConfig;
  private cameraState: CameraState;
  private isInitializedFlag: boolean = false;
  private snapshotCallbacks: Array<() => void> = [];

  // Lighting
  private ambientLight: THREE.AmbientLight | null = null;
  private directionalLight: THREE.DirectionalLight | null = null;

  /**
   * Create a new Three.js renderer
   *
   * @param container - DOM element to render into
   * @param config - Optional visualization configuration
   */
  constructor(container: HTMLElement, config?: Partial<VisualizationConfig>) {
    if (!container) {
      throw new Error('Container element is required for ThreeJsRenderer');
    }

    this.container = container;
    this.config = this.getDefaultConfig(config);
    this.cameraState = this.getDefaultCameraState();

    this.initialize();
  }

  /**
   * Initialize Three.js scene, camera, and renderer
   */
  private initialize(): void {
    // Create scene
    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(this.config.bgColor);

    // Create camera
    const aspect = this.container!.clientWidth / this.container!.clientHeight;
    this.camera = new THREE.PerspectiveCamera(75, aspect, 0.1, 1000);
    this.resetCamera();

    // Create renderer
    this.renderer = new THREE.WebGLRenderer({
      antialias: this.config.antialiasing,
      alpha: true,
    });
    this.updateRendererSize();

    // Add renderer to container
    this.container!.appendChild(this.renderer.domElement);

    // Initialize sub-renderers
    this.atomRenderer = new AtomRenderer(this.scene, this.config);
    this.bondRenderer = new BondRenderer(this.scene, this.config.bondRadius);
    this.cameraController = new CameraController(this.camera, this.cameraState);

    // Add lights
    this.setupLighting();

    // Add axes if enabled
    if (this.config.showAxes) {
      this.addAxesHelper();
    }

    // Start animation loop
    this.startAnimationLoop();

    this.isInitializedFlag = true;

    // Handle window resize
    window.addEventListener('resize', this.handleResize.bind(this));
  }

  /**
   * Setup lighting for the scene
   */
  private setupLighting(): void {
    if (!this.scene) {
      return;
    }

    this.ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
    this.scene.add(this.ambientLight);

    this.directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
    this.directionalLight.position.set(10, 10, 10);
    this.scene.add(this.directionalLight);
  }

  /**
   * Add axes helper to the scene
   */
  private addAxesHelper(): void {
    if (!this.scene) {
      return;
    }

    this.axesHelper = new THREE.AxesHelper(5);
    this.scene.add(this.axesHelper);
  }

  /**
   * Update renderer size based on container
   */
  private updateRendererSize(): void {
    if (!this.renderer || !this.container) {
      return;
    }

    const width = this.container.clientWidth;
    const height = this.container.clientHeight;

    this.renderer.setSize(width, height);

    if (this.camera) {
      this.camera.aspect = width / height;
      this.camera.updateProjectionMatrix();
    }
  }

  /**
   * Handle window resize event
   */
  private handleResize(): void {
    this.updateRendererSize();
  }

  /**
   * Start the animation loop
   */
  private startAnimationLoop(): void {
    const animate = () => {
      this.animationId = requestAnimationFrame(animate);

      if (this.renderer && this.scene && this.camera) {
        this.renderer.render(this.scene, this.camera);
      }
    };

    animate();
  }

  /**
   * Stop the animation loop
   */
  private stopAnimationLoop(): void {
    if (this.animationId !== null) {
      cancelAnimationFrame(this.animationId);
      this.animationId = null;
    }
  }

  /**
   * Get default visualization configuration
   */
  private getDefaultConfig(overrides?: Partial<VisualizationConfig>): VisualizationConfig {
    return {
      representationMode: 'ball-and-stick',
      showBonds: true,
      showUnitCell: false,
      showAxes: false,
      atomScale: 1.0,
      bondRadius: 0.1,
      bgColor: '#1a1a2e',
      antialiasing: true,
      ...overrides,
    };
  }

  /**
   * Get default camera state
   */
  private getDefaultCameraState(): CameraState {
    return {
      position: { x: 5, y: 5, z: 5 },
      target: { x: 0, y: 0, z: 0 },
      zoom: 1,
    };
  }

  /**
   * Check if renderer is initialized
   */
  public isInitialized(): boolean {
    return this.isInitializedFlag;
  }

  /**
   * Get the Three.js scene
   */
  public getScene(): THREE.Scene | null {
    return this.scene;
  }

  /**
   * Get the camera
   */
  public getCamera(): THREE.PerspectiveCamera | null {
    return this.camera;
  }

  /**
   * Get renderer size
   */
  public getSize(): { width: number; height: number } {
    if (!this.container) {
      return { width: 0, height: 0 };
    }

    return {
      width: this.container.clientWidth,
      height: this.container.clientHeight,
    };
  }

  /**
   * Render atoms to the scene
   *
   * @param atoms - Array of atoms to render
   * @returns Render result with success status and atom count
   */
  public renderAtoms(atoms: Atom[]): RenderResult {
    if (!this.scene || !this.camera || !this.atomRenderer) {
      return { success: false, atomCount: 0 };
    }

    // Store atoms
    this.atoms = atoms;

    const { atomCount, renderTime } = this.atomRenderer.renderAtoms(atoms);

    // Auto-center camera on structure
    this.cameraController?.centerOnAtoms(atoms);

    return {
      success: true,
      atomCount,
      renderTime,
    };
  }

  /**
   * Get current atom colors map
   */
  public getAtomColors(): Array<{ element: string; color: string }> {
    return this.atomRenderer?.getAtomColors(this.atoms) || [];
  }

  /**
   * Get current atom radii map
   */
  public getAtomRadii(): Record<string, number> {
    return this.atomRenderer?.getAtomRadii(this.atoms) || {};
  }

  /**
   * Get number of currently rendered atoms
   */
  public getAtomCount(): number {
    return this.atomRenderer?.getAtomCount() || 0;
  }

  /**
   * Detect and render bonds between atoms
   */
  public detectAndRenderBonds(): void {
    if (!this.bondRenderer) {
      return;
    }

    this.bonds = this.bondRenderer.detectAndRenderBonds(this.atoms, this.config.showBonds);
  }

  /**
   * Get number of currently rendered bonds
   */
  public getBondCount(): number {
    return this.bondRenderer?.getBondCount() || 0;
  }

  /**
   * Get bond threshold for two elements
   */
  public getBondThreshold(element1: string, element2: string): number {
    return BondRenderer.getBondThreshold(element1, element2);
  }

  /**
   * Render unit cell box
   *
   * @param lattice - 3x3 array of lattice vectors
   * @returns Render result
   */
  public renderUnitCell(lattice: number[][] | null): RenderResult {
    if (!this.scene || !lattice) {
      return { success: true, hasCell: false, atomCount: this.atoms.length };
    }

    // Clear existing unit cell
    this.clearUnitCell();

    // Create unit cell geometry
    const points: THREE.Vector3[] = [];

    const origin = new THREE.Vector3(0, 0, 0);
    const a = new THREE.Vector3(lattice[0][0], lattice[0][1], lattice[0][2]);
    const b = new THREE.Vector3(lattice[1][0], lattice[1][1], lattice[1][2]);
    const c = new THREE.Vector3(lattice[2][0], lattice[2][1], lattice[2][2]);

    const edges = [
      [origin, a],
      [origin, b],
      [origin, c],
      [a, a.clone().add(b)],
      [a, a.clone().add(c)],
      [b, b.clone().add(a)],
      [b, b.clone().add(c)],
      [c, c.clone().add(a)],
      [c, c.clone().add(b)],
      [a.clone().add(b), a.clone().add(b).add(c)],
      [a.clone().add(c), a.clone().add(b).add(c)],
      [b.clone().add(c), a.clone().add(b).add(c)],
    ];

    edges.forEach(([start, end]) => {
      points.push(start, end);
    });

    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    const material = new THREE.LineBasicMaterial({ color: 0x00ff00 });

    this.unitCellLines = new THREE.LineSegments(geometry, material);
    this.scene.add(this.unitCellLines);

    return { success: true, hasCell: true, atomCount: this.atoms.length };
  }

  /**
   * Clear unit cell from the scene
   */
  private clearUnitCell(): void {
    if (!this.scene || !this.unitCellLines) {
      return;
    }

    this.scene.remove(this.unitCellLines);
    this.unitCellLines.geometry.dispose();
    (this.unitCellLines.material as THREE.Material).dispose();
    this.unitCellLines = null;
  }

  /**
   * Rotate camera
   */
  public rotateCamera(deltaX: number, deltaY: number): void {
    this.cameraController?.rotate(deltaX, deltaY);
  }

  /**
   * Zoom camera
   */
  public zoomCamera(factor: number): void {
    this.cameraController?.zoom(factor);
  }

  /**
   * Pan camera
   */
  public panCamera(dx: number, dy: number, dz: number): void {
    this.cameraController?.pan(dx, dy, dz);
  }

  /**
   * Reset camera to default position
   */
  public resetCamera(): void {
    this.cameraState = this.getDefaultCameraState();

    if (this.camera) {
      this.camera.position.set(
        this.cameraState.position.x,
        this.cameraState.position.y,
        this.cameraState.position.z
      );
      this.camera.lookAt(
        this.cameraState.target.x,
        this.cameraState.target.y,
        this.cameraState.target.z
      );
    }

    if (this.cameraController) {
      this.cameraController.reset(this.cameraState);
    }
  }

  /**
   * Get camera zoom level
   */
  public getCameraZoom(): number {
    return this.cameraController?.getZoom() || 1;
  }

  /**
   * Update container size
   */
  public updateSize(width: number, height: number): void {
    if (!this.container) {
      return;
    }

    this.container.style.width = `${width}px`;
    this.container.style.height = `${height}px`;
    this.updateRendererSize();
  }

  /**
   * Set representation mode
   */
  public setRepresentationMode(mode: RepresentationMode): void {
    this.config.representationMode = mode;

    this.atomRenderer?.updateConfig(this.config);

    // Re-render atoms with new mode
    if (this.atoms.length > 0) {
      this.renderAtoms(this.atoms);
      if (this.config.showBonds) {
        this.detectAndRenderBonds();
      }
    }
  }

  /**
   * Get current representation mode
   */
  public getRepresentationMode(): RepresentationMode {
    return this.config.representationMode;
  }

  /**
   * Enable or disable lazy loading for large systems
   */
  public enableLazyLoading(enabled: boolean): void {
    // TODO: Implement lazy loading for large molecular systems
    // This will be part of Phase 5 optimization
  }

  /**
   * Export current scene as image
   */
  public exportImage(options: ExportOptions): { format: string; data: string } {
    if (!this.renderer) {
      throw new Error('Renderer not initialized');
    }

    const format = options.format || 'png';
    const mimeType = `image/${format}`;

    try {
      const dataUrl = this.renderer.domElement.toDataURL(mimeType, options.quality);
      return { format, data: dataUrl };
    } catch (error) {
      throw new Error(`Failed to export image: ${error}`);
    }
  }

  /**
   * Register callback for snapshot completion
   */
  public onSnapshotComplete(callback: () => void): void {
    this.snapshotCallbacks.push(callback);
  }

  /**
   * Take a snapshot of the current scene
   */
  public takeSnapshot(): void {
    if (this.renderer && this.scene && this.camera) {
      this.renderer.render(this.scene, this.camera);
    }

    this.snapshotCallbacks.forEach(cb => cb());
  }

  /**
   * Get current renderer state
   */
  public getState(): RendererState {
    return {
      atoms: this.atoms,
      bonds: this.bonds,
      camera: this.cameraController?.getState() || this.cameraState,
      config: this.config,
    };
  }

  /**
   * Dispose of all resources
   */
  public dispose(): void {
    if (!this.isInitializedFlag) {
      return;
    }

    // Stop animation
    this.stopAnimationLoop();

    // Clear sub-renderers
    this.atomRenderer?.clearAtoms();
    this.bondRenderer?.clearBonds();
    this.clearUnitCell();

    // Dispose lights
    if (this.ambientLight) {
      this.scene?.remove(this.ambientLight);
      this.ambientLight = null;
    }

    if (this.directionalLight) {
      this.scene?.remove(this.directionalLight);
      this.directionalLight = null;
    }

    if (this.axesHelper) {
      this.scene?.remove(this.axesHelper);
      this.axesHelper = null;
    }

    // Dispose renderer
    if (this.renderer) {
      if (this.container && this.renderer.domElement.parentNode) {
        this.container.removeChild(this.renderer.domElement);
      }
      this.renderer.dispose();
      this.renderer = null;
    }

    // Clear references
    this.scene = null;
    this.camera = null;
    this.container = null;
    this.atomRenderer = null;
    this.bondRenderer = null;
    this.cameraController = null;
    this.isInitializedFlag = false;

    // Remove event listeners
    window.removeEventListener('resize', this.handleResize.bind(this));
  }
}
