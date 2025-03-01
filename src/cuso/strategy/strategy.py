"""Abstract strategy class."""

from abc import abstractmethod
import logging
from typing import Any

class Strategy:
    """Base class to represent a Strategy.

    A Strategy represents an instance of a problem without an easily
    checked solution.
    """
    def __init__(self):
        self.logger = logging.getLogger("cuso." + self.__class__.__name__)

    @abstractmethod
    def run(self) -> Any:
        """Run the strategy."""