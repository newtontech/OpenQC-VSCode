#!/usr/bin/env python3
"""
Structure Writer CLI - Write structures from JSON to various formats
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from converter import ASEConverter


def main():
    parser = argparse.ArgumentParser(description='Write structures from JSON to various formats')
    parser.add_argument('command', choices=['write'], help='Command to execute')
    parser.add_argument('--format', '-f', required=True, help='Output format')
    parser.add_argument('--output', '-o', help='Output file path')
    parser.add_argument('--stdout', action='store_true', help='Write to stdout')
    parser.add_argument('--atoms-json', required=True, help='JSON string with atoms data')

    args = parser.parse_args()

    # Parse JSON input
    try:
        atoms_dict = json.loads(args.atoms_json)
    except json.JSONDecodeError as e:
        result = {'success': False, 'error': f'Invalid JSON input: {str(e)}'}
        print(json.dumps(result, indent=2))
        sys.exit(1)

    # Create converter and parse atoms from dict
    converter = ASEConverter()
    atoms_result = converter.atoms_from_dict(atoms_dict)

    if not atoms_result.success:
        result = {'success': False, 'error': atoms_result.error, 'warnings': atoms_result.warnings}
        print(json.dumps(result, indent=2))
        sys.exit(1)

    atoms = atoms_result.data

    # Write to stdout or file
    try:
        from ase.io import write

        if args.stdout:
            write('-', atoms, format=args.format)
            result = {'success': True, 'output_file': 'stdout', 'metadata': {'output_format': args.format, 'natoms': len(atoms)}}
        else:
            if not args.output:
                result = {'success': False, 'error': 'No output path specified'}
                print(json.dumps(result, indent=2))
                sys.exit(1)

            write(args.output, atoms, format=args.format)
            result = {'success': True, 'data': args.output, 'metadata': {'output_file': args.output, 'output_format': args.format, 'natoms': len(atoms)}}

        print(json.dumps(result, indent=2))

    except Exception as e:
        result = {'success': False, 'error': f'Failed to write structure: {str(e)}'}
        print(json.dumps(result, indent=2))
        sys.exit(1)


if __name__ == '__main__':
    main()
