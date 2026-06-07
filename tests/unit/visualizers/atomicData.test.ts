import {
  ELEMENT_SYMBOLS_BY_ATOMIC_NUMBER,
  UNKNOWN_ELEMENT_SYMBOL,
  getElementSymbol,
} from '../../../src/visualizers/atomicData';

describe('atomicData', () => {
  it('contains symbols for all currently named elements', () => {
    expect(ELEMENT_SYMBOLS_BY_ATOMIC_NUMBER).toHaveLength(119);
    expect(getElementSymbol(1)).toBe('H');
    expect(getElementSymbol(54)).toBe('Xe');
    expect(getElementSymbol(78)).toBe('Pt');
    expect(getElementSymbol(92)).toBe('U');
    expect(getElementSymbol(118)).toBe('Og');
  });

  it('falls back for invalid or unsupported atomic numbers', () => {
    expect(getElementSymbol(0)).toBe(UNKNOWN_ELEMENT_SYMBOL);
    expect(getElementSymbol(119)).toBe(UNKNOWN_ELEMENT_SYMBOL);
    expect(getElementSymbol(Number.NaN)).toBe(UNKNOWN_ELEMENT_SYMBOL);
    expect(getElementSymbol(1.5)).toBe(UNKNOWN_ELEMENT_SYMBOL);
  });
});
