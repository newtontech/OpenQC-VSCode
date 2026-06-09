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

/** Full backend check result. */
export interface BackendCheckResult {
  success: boolean;
  python: PythonInfo;
  packages: Record<string, PackageStatus>;
  externalTools: Record<string, ExternalToolStatus>;
}

/** Bridge error from Python subprocess. */
export interface BridgeError {
  message: string;
  [key: string]: unknown;
}
