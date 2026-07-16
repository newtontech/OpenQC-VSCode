/**
 * Types and interfaces for the Scientific Python backend status.
 *
 * @module python/PythonBackendStatus
 */

/** Information about the Python interpreter. */
export interface PythonInfo {
  executable: string;
  version: string;
  platform?: string;
}

/** Status of a single Python package. */
export interface PackageStatus {
  available: boolean;
  version?: string;
  installHint?: string;
}

/** Status of an external tool. */
export interface ExternalToolStatus {
  available: boolean;
  path: string | null;
}

export type BackendHealthStatus = 'installed' | 'degraded' | 'missing';
export type BackendCapabilityStatus = 'available' | 'degraded' | 'missing';

export interface BackendCapability {
  label?: string;
  status: BackendCapabilityStatus;
  detail: string;
  requires?: string[];
}

/** Full backend check result. */
export interface BackendCheckResult {
  success: boolean;
  status?: BackendHealthStatus;
  statusDetail?: string;
  python: PythonInfo;
  packages: Record<string, PackageStatus>;
  externalTools: Record<string, ExternalToolStatus>;
  capabilities?: Record<string, BackendCapability>;
  missingPackages?: string[];
}

/** Bridge error from Python subprocess. */
export interface BridgeError {
  message: string;
  [key: string]: unknown;
}
