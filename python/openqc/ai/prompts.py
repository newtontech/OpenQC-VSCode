"""
Prompt Templates for Quantum Chemistry AI

Provides context-aware prompt templates for various quantum chemistry
codes and AI tasks.
"""

from typing import Dict, Optional, Any


class PromptTemplates:
    """Quantum chemistry prompt templates."""

    # Software-specific context
    SOFTWARE_CONTEXT = {
        'vasp': """
VASP (Vienna Ab initio Simulation Package) is a plane-wave DFT code.
Key concepts: INCAR (parameters), POSCAR (structure), KPOINTS (k-grid), POTCAR (pseudopotentials).
Important parameters: ENCUT (cutoff), ISPIN (spin), ISMEAR (smearing), SIGMA (smearing width).
""",
        'cp2k': """
CP2K is a hybrid Gaussian and plane-wave DFT code.
Key sections: GLOBAL, FORCE_EVAL, MOTION, CELL, COORD, KIND.
Important parameters: CUTOFF, REL_CUTOFF, BASIS_SET_FILE_NAME, POTENTIAL_FILE_NAME.
""",
        'qe': """
Quantum ESPRESSO is a plane-wave DFT code.
Key namelists: CONTROL, SYSTEM, ELECTRONS, IONS, CELL.
Important parameters: ecutwfc, ecutrho, conv_thr, mixing_beta.
""",
        'gaussian': """
Gaussian is a quantum chemistry code using Gaussian basis sets.
Key sections: route card, title, molecule specification.
Important keywords: Opt, Freq, MP2, B3LYP, 6-31G(d).
""",
        'orca': """
ORCA is a quantum chemistry code with emphasis on spectroscopy.
Key sections: method, basis, job type.
Important keywords: B3LYP, def2-TZVP, Opt, Freq, TightSCF.
""",
        'nwchem': """
NWChem is a computational chemistry code.
Key sections: geometry, basis, scf, dft, task.
Important keywords: basis, dft, task, memory.
""",
        'gamess': """
GAMESS is a quantum chemistry code.
Key sections: $CONTRL, $BASIS, $DATA, $SYSTEM.
Important parameters: RUNTYP, DFTTYP, GBASIS.
""",
    }

    @classmethod
    def get_optimize_prompt(
        cls,
        input_content: str,
        software: str,
        context: Optional[Dict[str, Any]] = None
    ) -> str:
        """Generate prompt for input optimization."""
        software_context = cls.SOFTWARE_CONTEXT.get(software, '')
        
        prompt = f"""You are an expert in quantum chemistry computational methods.

{software_context}

Please analyze the following {software.upper()} input file and suggest optimizations.
Focus on:
1. Accuracy vs efficiency trade-offs
2. Convergence settings
3. Parameter consistency
4. Missing important parameters
5. Best practices

Input file:
```
{input_content}
```

Provide your response in the following JSON format:
{{
    "suggestions": [
        {{
            "type": "optimization",
            "parameter": "parameter_name",
            "currentValue": "current_value",
            "suggestedValue": "suggested_value",
            "message": "brief description",
            "explanation": "detailed explanation"
        }}
    ],
    "content": "overall assessment"
}}
"""
        return prompt

    @classmethod
    def get_generate_prompt(
        cls,
        description: str,
        software: str,
        context: Optional[Dict[str, Any]] = None
    ) -> str:
        """Generate prompt for input generation."""
        software_context = cls.SOFTWARE_CONTEXT.get(software, '')
        
        structure_info = ""
        if context and 'structure' in context:
            structure_info = f"""
Structure information:
```
{context['structure']}
```
"""
        
        prompt = f"""You are an expert in quantum chemistry computational methods.

{software_context}

Please generate a {software.upper()} input file for the following calculation:

Description: {description}

{structure_info}

Generate a complete, ready-to-run input file with appropriate parameters.
Provide only the input file content without any additional explanation.
"""
        return prompt

    @classmethod
    def get_explain_prompt(
        cls,
        input_content: str,
        software: str
    ) -> str:
        """Generate prompt for parameter explanation."""
        software_context = cls.SOFTWARE_CONTEXT.get(software, '')
        
        prompt = f"""You are an expert in quantum chemistry computational methods.

{software_context}

Please explain the following {software.upper()} input file parameters in plain English:

Input file:
```
{input_content}
```

For each parameter/section:
1. Explain what it does
2. Explain why it's set this way
3. Mention typical values and their implications
4. Note any potential issues or considerations

Provide a clear, educational explanation suitable for someone learning computational chemistry.
"""
        return prompt

    @classmethod
    def get_debug_prompt(
        cls,
        input_content: str,
        output_content: str,
        software: str
    ) -> str:
        """Generate prompt for debugging."""
        software_context = cls.SOFTWARE_CONTEXT.get(software, '')
        
        prompt = f"""You are an expert in quantum chemistry computational methods and debugging.

{software_context}

A calculation has failed or produced unexpected results. Please analyze the input and output files:

Input file:
```
{input_content}
```

Output/log file:
```
{output_content}
```

Please:
1. Identify any errors or warnings
2. Suggest the likely cause of the problem
3. Recommend specific fixes
4. Provide guidance on preventing similar issues

Be specific and actionable in your recommendations.
"""
        return prompt

    @classmethod
    def get_software_list(cls) -> list:
        """Get list of supported software."""
        return list(cls.SOFTWARE_CONTEXT.keys())
