"""
   Copyright 2024 Manas Mahale

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

__version__ = "0.1.0"
__author__ = "Manas Mahale"
__email__ = "manas.mahale@bcp.edu.in"
__license__ = "Apache License 2.0"
__url__ = "https://github.com/Manas02/fbdd"
__keywords__ = ["AI", "Machine Learning", "BERT", "Cheminformatics"]

from .molecule import Molecule
from .error import InvalidSMILESError, InvalidSMARTSError
from .fragment import Fragment
