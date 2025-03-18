# Standard libraries
from dataclasses import dataclass
from typing import ClassVar


@dataclass
class Settings:
    """
    This data-class is used to define and store plot settings which are not provided as arguments.
    This includes things like the figure dimensions, colors, symbols font sizes etc.
    """

    symbols: ClassVar[list] = ["o", "v", "s", "^", "d", "p"]
    colors: ClassVar[list] = ["blue", "orange", "limegreen", "red"]
    figsize: ClassVar[tuple] = (3.5, 2.3)
    axes_size: int = 9
    tick_size: int = 9
