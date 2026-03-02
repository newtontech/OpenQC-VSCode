/**
 * Three.js-based 3D molecular structure renderer
 *
 * Phase 2: 3D Visualization - Core Renderer
 *
 * This module provides Three.js integration for rendering molecular structures
 * from quantum chemistry input files (VASP POSCAR, Gaussian, ORCA, etc.)
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
  ELEMENT_COLORS,
  COVALENT_RADII,
  VDW_RADII,
} from './types';

/**
 * Three.js-based molecular structure renderer
 */
export class ThreeJsRenderer {
  private scene: THREE.Scene | null = null;
  private camera: THREE.PerspectiveCamera | null = null;
  private renderer: THREE.WebGLRenderer | null = null;
  private container: HTMLElement | null = null;
  private animationId: number | null = null;

  // 3D objects storage
  private atomMeshes: Map<number, THREE.Mesh> = new Map();
  private bondMeshes: Map<string, THREE.Mesh> = new Map();
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
    if (!this.scene) return;

    // Ambient light for overall illumination
    this.ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
    this.scene.add(this.ambientLight);

    // Directional light for shadows and depth
    this.directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
    this.directionalLight.position.set(10, 10, 10);
    this.scene.add(this.directionalLight);
  }

  /**
   * Add axes helper to the scene
   */
  private addAxesHelper(): void {
    if (!this.scene) return;

    this.axesHelper = new THREE.AxesHelper(5);
    this.scene.add(this.axesHelper);
  }

  /**
   * Update renderer size based on container
   */
  private updateRendererSize(): void {
    if (!this.renderer || !this.container) return;

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
  private getDefaultConfig(
    overrides?: Partial<VisualizationConfig>
  ): VisualizationConfig {
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
    if (!this.scene || !this.camera) {
      return { success: false, atomCount: 0 };
    }

    const startTime = performance.now();

    // Clear existing atoms
    this.clearAtoms();

    // Store atoms
    this.atoms = atoms;

    // Create geometry and material for each atom
    atoms.forEach((atom, index) => {
      const mesh = this.createAtomMesh(atom);
      if (mesh) {
        this.atomMeshes.set(index, mesh);
        this.scene!.add(mesh);
      }
    });

    const renderTime = performance.now() - startTime;

    // Auto-center camera on structure
    this.centerCameraOnStructure();

    return {
      success: true,
      atomCount: atoms.length,
      renderTime,
    };
  }

  /**
   * Create a mesh for a single atom
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
   * Get radius for an atom based on representation mode
   */
  private getAtomRadius(element: string): number {
    const baseRadius =
      this.config.representationMode === 'space-filling'
        ? VDW_RADII[element] || VDW_RADII.default
        : COVALENT_RADII[element] || COVALENT_RADII.default;

    return baseRadius * this.config.atomScale;
  }

  /**
   * Get color for an element
   */
  private getAtomColor(element: string): string {
    return ELEMENT_COLORS[element] || ELEMENT_COLORS.default;
  }

  /**
   * Get current atom colors map
   */
  public getAtomColors(): Array<{ element: string; color: string }> {
    const uniqueElements = [...new Set(this.atoms.map((a) => a.element))];
    return uniqueElements.map((element) => ({
      element,
      color: this.getAtomColor(element),
    }));
  }

  /**
   * Get current atom radii map
   */
  public getAtomRadii(): Record<string, number> {
    const uniqueElements = [...new Set(this.atoms.map((a) => a.element))];
    const radii: Record<string, number> = {};
    uniqueElements.forEach((element) => {
      radii[element] = this.getAtomRadius(element);
    });
    return radii;
  }

  /**
   * Get number of currently rendered atoms
   */
  public getAtomCount(): number {
    return this.atomMeshes.size;
  }

  /**
   * Clear all atom meshes from the scene
   */
  private clearAtoms(): void {
    if (!this.scene) return;

    this.atomMeshes.forEach((mesh) => {
      this.scene!.remove(mesh);
      mesh.geometry.dispose();
      if (Array.isArray(mesh.material)) {
        mesh.material.forEach((m) => m.dispose());
      } else {
        mesh.material.dispose();
      }
    });

    this.atomMeshes.clear();
  }

  /**
   * Detect and render bonds between atoms
   */
  public detectAndRenderBonds(): void {
    if (!this.config.showBonds) return;

    // Clear existing bonds
    this.clearBonds();

    // Detect bonds
    this.bonds = this.detectBonds();

    // Render bonds
    this.bonds.forEach((bond, index) => {
      const mesh = this.createBondMesh(bond);
      if (mesh) {
        this.bondMeshes.set(`${bond.atomIndex1}-${bond.atomIndex2}`, mesh);
        this.scene!.add(mesh);
      }
    });
  }

  /**
   * Detect bonds between atoms based on distance
   */
  private detectBonds(): Bond[] {
    const bonds: Bond[] = [];
    const tolerance = 1.2; // 20% tolerance for bond detection

    for (let i = 0; i < this.atoms.length; i++) {
      for (let j = i + 1; j < this.atoms.length; j++) {
        const atom1 = this.atoms[i];
        const atom2 = this.atoms[j];

        const distance = this.calculateDistance(atom1, atom2);
        const maxDistance =
          (COVALENT_RADII[atom1.element] || COVALENT_RADII.default) +
          (COVALENT_RADII[atom2.element] || COVALENT_RADII.default);

        if (distance <= maxDistance * tolerance) {
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
   * Calculate distance between two atoms
   */
  private calculateDistance(atom1: Atom, atom2: Atom): number {
    const dx = atom2.x - atom1.x;
    const dy = atom2.y - atom1.y;
    const dz = atom2.z - atom1.z;

    return Math.sqrt(dx * dx + dy * dy + dz * dz);
  }

  /**
   * Create a mesh for a single bond
   */
  private createBondMesh(bond: Bond): THREE.Mesh | null {
    const atom1 = this.atoms[bond.atomIndex1];
    const atom2 = this.atoms[bond.atomIndex2];

    if (!atom1 || !atom2) return null;

    // Create cylinder for bond
    const direction = new THREE.Vector3(
      atom2.x - atom1.x,
      atom2.y - atom1.y,
      atom2.z - atom1.z
    );
    const length = direction.length();

    const geometry = new THREE.CylinderGeometry(
      this.config.bondRadius,
      this.config.bondRadius,
      length,
      8
    );

    const material = new THREE.MeshPhongMaterial({
      color: 0x888888,
    });

    const mesh = new THREE.Mesh(geometry, material);

    // Position and orient the cylinder
    const midpoint = new THREE.Vector3(
      (atom1.x + atom2.x) / 2,
      (atom1.y + atom2.y) / 2,
      (atom1.z + atom2.z) / 2
    );
    mesh.position.copy(midpoint);
    mesh.quaternion.setFromUnitVectors(
      new THREE.Vector3(0, 1, 0),
      direction.normalize()
    );

    return mesh;
  }

  /**
   * Get number of currently rendered bonds
   */
  public getBondCount(): number {
    return this.bondMeshes.size;
  }

  /**
   * Clear all bond meshes from the scene
   */
  private clearBonds(): void {
    if (!this.scene) return;

    this.bondMeshes.forEach((mesh) => {
      this.scene!.remove(mesh);
      mesh.geometry.dispose();
      if (Array.isArray(mesh.material)) {
        mesh.material.forEach((m) => m.dispose());
      } else {
        mesh.material.dispose();
      }
    });

    this.bondMeshes.clear();
  }

  /**
   * Get bond threshold for two elements
   */
  public getBondThreshold(element1: string, element2: string): number {
    return (
      (COVALENT_RADII[element1] || COVALENT_RADII.default) +
      (COVALENT_RADII[element2] || COVALENT_RADII.default)
    );
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

    // Origin and corners
    const origin = new THREE.Vector3(0, 0, 0);
    const a = new THREE.Vector3(lattice[0][0], lattice[0][1], lattice[0][2]);
    const b = new THREE.Vector3(lattice[1][0], lattice[1][1], lattice[1][2]);
    const c = new THREE.Vector3(lattice[2][0], lattice[2][1], lattice[2][2]);

    // Edges
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
    if (!this.scene || !this.unitCellLines) return;

    this.scene.remove(this.unitCellLines);
    this.unitCellLines.geometry.dispose();
    (this.unitCellLines.material as THREE.Material).dispose();
    this.unitCellLines = null;
  }

  /**
   * Center camera on the molecular structure
   */
  private centerCameraOnStructure(): void {
    if (!this.camera || this.atoms.length === 0) return;

    // Calculate bounding box
    const box = new THREE.Box3();
    this.atoms.forEach((atom) => {
      box.expandByPoint(new THREE.Vector3(atom.x, atom.y, atom.z));
    });

    // Center the box
    const center = box.getCenter(new THREE.Vector3());

    // Update camera target
    this.cameraState.target = { x: center.x, y: center.y, z: center.z };

    // Position camera
    const size = box.getSize(new THREE.Vector3());
    const maxDim = Math.max(size.x, size.y, size.z);
    const distance = maxDim * 2;

    this.camera.position.set(
      center.x + distance,
      center.y + distance,
      center.z + distance
    );
    this.camera.lookAt(center);
  }

  /**
   * Rotate camera
   */
  public rotateCamera(deltaX: number, deltaY: number): void {
    if (!this.camera) return;

    // Simple orbital rotation around target
    const theta = (deltaX * Math.PI) / 180;
    const phi = (deltaY * Math.PI) / 180;

    // Get current position relative to target
    const offset = new THREE.Vector3(
      this.camera.position.x - this.cameraState.target.x,
      this.camera.position.y - this.cameraState.target.y,
      this.camera.position.z - this.cameraState.target.z
    );

    // Rotate around Y axis
    offset.applyAxisAngle(new THREE.Vector3(0, 1, 0), theta);

    // Rotate around X axis (limited)
    const xAxis = new THREE.Vector3(1, 0, 0).applyAxisAngle(
      new THREE.Vector3(0, 1, 0),
      theta
    );
    offset.applyAxisAngle(xAxis, phi);

    // Update camera position
    this.camera.position.set(
      this.cameraState.target.x + offset.x,
      this.cameraState.target.y + offset.y,
      this.cameraState.target.z + offset.z
    );

    this.camera.lookAt(
      this.cameraState.target.x,
      this.cameraState.target.y,
      this.cameraState.target.z
    );
  }

  /**
   * Zoom camera
   */
  public zoomCamera(factor: number): void {
    if (!this.camera) return;

    this.cameraState.zoom *= factor;

    const offset = new THREE.Vector3(
      this.camera.position.x - this.cameraState.target.x,
      this.camera.position.y - this.cameraState.target.y,
      this.camera.position.z - this.cameraState.target.z
    );

    offset.multiplyScalar(factor);

    this.camera.position.set(
      this.cameraState.target.x + offset.x,
      this.cameraState.target.y + offset.y,
      this.cameraState.target.z + offset.z
    );
  }

  /**
   * Pan camera
   */
  public panCamera(dx: number, dy: number, dz: number): void {
    if (!this.camera) return;

    this.camera.position.x += dx;
    this.camera.position.y += dy;
    this.camera.position.z += dz;

    this.cameraState.target.x += dx;
    this.cameraState.target.y += dy;
    this.cameraState.target.z += dz;
  }

  /**
   * Reset camera to default position
   */
  public resetCamera(): void {
    if (!this.camera) return;

    this.cameraState = this.getDefaultCameraState();

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

  /**
   * Get camera zoom level
   */
  public getCameraZoom(): number {
    return this.cameraState.zoom;
  }

  /**
   * Update container size
   */
  public updateSize(width: number, height: number): void {
    if (!this.container) return;

    this.container.style.width = `${width}px`;
    this.container.style.height = `${height}px`;
    this.updateRendererSize();
  }

  /**
   * Set representation mode
   */
  public setRepresentationMode(mode: RepresentationMode): void {
    this.config.representationMode = mode;

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
    // Trigger render and notify callbacks
    if (this.renderer && this.scene && this.camera) {
      this.renderer.render(this.scene, this.camera);
    }

    this.snapshotCallbacks.forEach((cb) => cb());
  }

  /**
   * Get current renderer state
   */
  public getState(): RendererState {
    return {
      atoms: this.atoms,
      bonds: this.bonds,
      camera: this.cameraState,
      config: this.config,
    };
  }

  /**
   * Dispose of all resources
   */
  public dispose(): void {
    if (!this.isInitializedFlag) return;

    // Stop animation
    this.stopAnimationLoop();

    // Clear all objects
    this.clearAtoms();
    this.clearBonds();
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
    this.isInitializedFlag = false;

    // Remove event listeners
    window.removeEventListener('resize', this.handleResize.bind(this));
  }
}
