import * as path from 'path';
import * as vscode from 'vscode';
import type { QuantumChemistrySoftware } from '../managers/FileTypeDetector';
import { CacheKeyGenerator, CacheManager, type LRUCache } from '../performance/cacheManager';
import { IncrementalParser } from '../performance/incrementalParser';
import { getWorkerManager, type WorkerManager } from '../performance/workerManager';
import { WorkerMessageType } from '../performance/computeWorker';
import type { ASEAtoms } from '../ase/ASEConverter';
import {
  OPENQC_STRUCTURE_SCHEMA_VERSION,
  type OpenQCStructure,
} from '../structures/OpenQCStructure';
import {
  createOpenQCStructure,
  poscarToOpenQCStructure,
  molecularStructureToOpenQCStructure,
} from '../structures/converters';
import { StructureConverter } from './StructureConverter';
import { Molecule3D } from './Molecule3D';

export type ViewerParseSource = 'cache' | 'worker' | 'incremental';

export interface ViewerParseResult {
  structure: OpenQCStructure;
  source: ViewerParseSource;
  incremental: boolean;
}

/** Production path that connects worker parsing, LRU caching, and incremental parsing to 3Dmol. */
export class ViewerStructurePipeline {
  private readonly parsers = new Map<string, IncrementalParser<OpenQCStructure>>();
  private readonly structureCache: LRUCache<OpenQCStructure>;

  constructor(
    context: vscode.ExtensionContext,
    private readonly workerManager: WorkerManager = getWorkerManager()
  ) {
    this.structureCache = CacheManager.getInstance(context).getStructureCache();
  }

  async parse(
    document: vscode.TextDocument,
    software: QuantumChemistrySoftware
  ): Promise<ViewerParseResult> {
    const filePath = document.uri.fsPath || document.fileName || document.uri.toString();
    const cacheKey = CacheKeyGenerator.forFile(filePath, `viewer:${document.version}`);
    const cached = this.structureCache.get(cacheKey);
    if (cached) {
      return { structure: cached, source: 'cache', incremental: false };
    }

    const workerStructure = await this.tryWorkerParse(document, software);
    if (workerStructure) {
      this.structureCache.set(cacheKey, workerStructure);
      return { structure: workerStructure, source: 'worker', incremental: false };
    }

    const parserKey = `${document.uri.toString()}::${software}`;
    let parser = this.parsers.get(parserKey);
    if (!parser) {
      parser = new IncrementalParser(content =>
        parseViewerStructureContent(content, document.fileName, software)
      );
      this.parsers.set(parserKey, parser);
    }

    const parsed = await parser.parse(document);
    this.structureCache.set(cacheKey, parsed.ast);
    const incremental =
      parsed.changes.added.length > 0 ||
      parsed.changes.removed.length > 0 ||
      parsed.changes.modified.length > 0;
    return { structure: parsed.ast, source: 'incremental', incremental };
  }

  clear(): void {
    this.parsers.forEach(parser => parser.clearCache());
    this.parsers.clear();
  }

  private async tryWorkerParse(
    document: vscode.TextDocument,
    software: QuantumChemistrySoftware
  ): Promise<OpenQCStructure | undefined> {
    const basename = path.basename(document.fileName).toUpperCase();
    if (software !== 'VASP' || (basename !== 'POSCAR' && basename !== 'CONTCAR')) {
      return undefined;
    }

    try {
      const task = await this.workerManager.submitTask(
        WorkerMessageType.PARSE_STRUCTURE,
        { content: document.getText(), format: 'poscar' },
        'high'
      );
      const completed = await this.workerManager.waitForTask(task.id, 10_000);
      return aseAtomsToOpenQCStructure(completed.result as ASEAtoms, document.fileName, software);
    } catch {
      return undefined;
    }
  }
}

export function parseViewerStructureContent(
  content: string,
  fileName: string,
  software: QuantumChemistrySoftware
): OpenQCStructure {
  if (software === 'VASP') {
    const basename = path.basename(fileName).toUpperCase();
    if (basename === 'POSCAR' || basename === 'CONTCAR') {
      try {
        const structure = poscarToOpenQCStructure(content, fileName);
        if (structure.atoms.length > 0) {
          return structure;
        }
      } catch {
        // Continue through shared native fallbacks.
      }
    }
  }

  try {
    const molecular = new StructureConverter().autoConvert(content, fileName);
    if (molecular.atoms.length > 0) {
      return molecularStructureToOpenQCStructure(molecular, {
        sourceSoftware: software,
        sourceParser: 'native',
      });
    }
  } catch {
    // Continue to legacy format parsing.
  }

  const atoms = new Molecule3D().parseAtoms(content, software);
  if (atoms.length === 0) {
    throw new Error('No molecular structure found in file');
  }
  return createOpenQCStructure(atoms, { name: fileName, sourceSoftware: software });
}

function aseAtomsToOpenQCStructure(
  atoms: ASEAtoms,
  fileName: string,
  software: QuantumChemistrySoftware
): OpenQCStructure {
  const cell = atoms.cell
    ? {
        a: tuple3(atoms.cell[0]),
        b: tuple3(atoms.cell[1]),
        c: tuple3(atoms.cell[2]),
        pbc: tupleBoolean3(atoms.pbc),
        coordinateMode: 'cartesian' as const,
      }
    : undefined;

  return {
    schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
    kind: cell ? 'periodic' : 'molecule',
    name: fileName,
    atoms: atoms.chemical_symbols.map((element, index) => ({
      element,
      x: atoms.positions[index][0],
      y: atoms.positions[index][1],
      z: atoms.positions[index][2],
    })),
    cell,
    metadata: {
      source: { filename: fileName, software, parser: 'compute-worker-native' },
      provenance: { createdAt: new Date().toISOString(), warnings: [] },
    },
  };
}

function tuple3(values: number[]): [number, number, number] {
  return [values[0], values[1], values[2]];
}

function tupleBoolean3(values: boolean[]): [boolean, boolean, boolean] {
  return [Boolean(values[0]), Boolean(values[1]), Boolean(values[2])];
}
