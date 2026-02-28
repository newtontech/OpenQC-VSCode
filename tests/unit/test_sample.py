"""
Sample Python Unit Test

This demonstrates the basic structure for pytest unit tests.
"""

import pytest


def test_basic_assertion():
    """Test basic assertion."""
    assert True


def test_async_function():
    """Test async functionality."""
    result = 42
    assert result == 42


def test_error_handling():
    """Test that errors are raised correctly."""
    with pytest.raises(ValueError):
        raise ValueError("Test error")


class TestClassExample:
    """Example of using test classes."""
    
    def test_method(self):
        """Test method in class."""
        assert 1 + 1 == 2
    
    @pytest.mark.parametrize("input,expected", [
        (1, 2),
        (2, 4),
        (3, 6),
    ])
    def test_parametrized(self, input, expected):
        """Test with parameters."""
        assert input * 2 == expected
