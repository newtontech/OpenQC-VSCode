"""
Response Parser for AI Outputs

Parses LLM responses into structured data.
"""

import re
import json
from typing import Dict, List, Any, Optional


class ResponseParser:
    """Parse AI responses into structured formats."""

    def parse_optimization(self, content: str) -> Dict[str, Any]:
        """Parse optimization response."""
        result = {
            'suggestions': [],
            'content': content
        }
        
        # Try to extract JSON
        json_match = re.search(r'```json\s*(.*?)\s*```', content, re.DOTALL)
        if json_match:
            try:
                data = json.loads(json_match.group(1))
                result['suggestions'] = data.get('suggestions', [])
                result['content'] = data.get('content', content)
                return result
            except json.JSONDecodeError:
                pass
        
        # Try to parse inline JSON
        try:
            # Look for JSON object
            json_start = content.find('{')
            json_end = content.rfind('}')
            if json_start != -1 and json_end != -1:
                data = json.loads(content[json_start:json_end+1])
                result['suggestions'] = data.get('suggestions', [])
        except json.JSONDecodeError:
            pass
        
        return result

    def parse_generated_input(self, content: str) -> str:
        """Extract generated input file from response."""
        # Look for code blocks
        code_match = re.search(r'```\w*\s*(.*?)\s*```', content, re.DOTALL)
        if code_match:
            return code_match.group(1).strip()
        
        # If no code block, return content after a marker
        markers = [
            'Input file:',
            'Generated input:',
            'Here is the input file:',
        ]
        for marker in markers:
            idx = content.find(marker)
            if idx != -1:
                return content[idx + len(marker):].strip()
        
        # Return whole content
        return content.strip()

    def parse_explanation(self, content: str) -> List[Dict[str, str]]:
        """Parse explanation into sections."""
        sections = []
        
        # Split by headers
        lines = content.split('\n')
        current_section = None
        current_content = []
        
        for line in lines:
            # Check for header
            if line.startswith('#') or line.startswith('**'):
                if current_section:
                    sections.append({
                        'title': current_section,
                        'content': '\n'.join(current_content).strip()
                    })
                current_section = line.lstrip('#*').strip()
                current_content = []
            else:
                current_content.append(line)
        
        # Add last section
        if current_section:
            sections.append({
                'title': current_section,
                'content': '\n'.join(current_content).strip()
            })
        
        return sections

    def parse_debug(self, content: str) -> Dict[str, Any]:
        """Parse debug response."""
        result = {
            'errors': [],
            'warnings': [],
            'suggestions': [],
            'summary': content
        }
        
        # Extract errors
        error_pattern = r'(?:error|ERROR|Error)[s]?[:\s]*(.+?)(?=\n\n|\n(?:warning|WARNING)|$)'
        errors = re.findall(error_pattern, content, re.DOTALL | re.IGNORECASE)
        result['errors'] = [e.strip() for e in errors if e.strip()]
        
        # Extract warnings
        warning_pattern = r'(?:warning|WARNING|Warning)[s]?[:\s]*(.+?)(?=\n\n|\n(?:error|ERROR)|$)'
        warnings = re.findall(warning_pattern, content, re.DOTALL | re.IGNORECASE)
        result['warnings'] = [w.strip() for w in warnings if w.strip()]
        
        # Extract suggestions
        suggestion_patterns = [
            r'(?:suggest|recommend|fix)[ei]?[sd]?[:\s]*(.+?)(?=\n\n|$)',
            r'(?:try|consider)[\s:]+(.+?)(?=\n\n|$)',
        ]
        for pattern in suggestion_patterns:
            suggestions = re.findall(pattern, content, re.DOTALL | re.IGNORECASE)
            result['suggestions'].extend([s.strip() for s in suggestions if s.strip()])
        
        return result
