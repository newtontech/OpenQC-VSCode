export interface Atom {
  elem: string;
  x: number;
  y: number;
  z: number;
}

export interface UnitCell {
  a: number;
  b: number;
  c: number;
  alpha: number;
  beta: number;
  gamma: number;
}

export interface NGLStructureData {
  atoms: Atom[];
  unitCell?: UnitCell;
}

export class StructureConverter {
  /**
   * Convert atoms array to XYZ format string
   * @param atoms - Array of atoms with element symbol and coordinates
   * @param comment - Comment line for the XYZ file (usually molecule name)
   * @returns XYZ format string
   */
  static atomsToXYZ(atoms: Atom[], comment: string = 'molecule'): string {
    const lines: string[] = [];

    // First line: number of atoms
    lines.push(String(atoms.length));

    // Second line: comment
    lines.push(comment);

    // Atom lines: element x y z (6 decimal places)
    for (const atom of atoms) {
      const x = atom.x.toFixed(6);
      const y = atom.y.toFixed(6);
      const z = atom.z.toFixed(6);
      lines.push(`${atom.elem} ${x} ${y} ${z}`);
    }

    return lines.join('\n');
  }

  /**
   * Convert atoms array to NGL-compatible JSON format
   * @param atoms - Array of atoms with element symbol and coordinates
   * @param unitCell - Optional unit cell parameters for periodic systems
   * @returns JSON string
   */
  static atomsToJSON(atoms: Atom[], unitCell?: UnitCell): string {
    const data: NGLStructureData = { atoms };

    if (unitCell) {
      data.unitCell = unitCell;
    }

    return JSON.stringify(data);
  }
}
