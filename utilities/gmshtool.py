#!python

"""
To record a macro to some register, for example 0,
press 'q0' while in command mode and
start performing actions over the text.
To stop recording go to the command mode
and press 'q'
To reproduce macro use the '@0' command, or
just '@@' to reproduce the last macro
"""

import re

class GmshSection(object):
  """ 
  Class defining 2 items:
  1) a section name
  2) an array of lines in section
  """
  def __init__(self,section_string = ""):
    super(GmshSection,self).__init__()
    self.lines = []
    self.name = ""
    if len(section_string):
      self.load_section(section_string)

  def load_section(self,section_string):
    section = section_string.split('\n')
    self.name = section[0][1:]
    self.lines = section[1:-1]


class GmshSectionParser(object):
  """
  Base class for section parsers
  """
  section_name = ""
  def __init__(self):
    super(GmshSectionParser,self).__init__()

  def parse_section(self,section):
    assert(section.name == self.section_name)
    self.do_parse(section)

  def do_parse(self,section):
    self.values = section.lines

class GmshMeshFormat(GmshSectionParser):
  section_name = 'MeshFormat' 

class GmshNodes(GmshSectionParser):
  section_name = 'Nodes'
  def do_parse(self,section):
    self.nodes = []
    count = int(section.lines[0])
    self.nodes = map(lambda l:[float(i) for i in l.split()[1:]],
      section.lines[1:])

class GmshElements(GmshSectionParser):
  section_name = 'Elements'
  def do_parse(self,section):
    self.elements = []
    count = int(section.lines[0])
    self.elements = map(lambda l:[int(i)-1 for i in l.split()[6:]],
      section.lines[1:])

class GmshNodeData(GmshSectionParser):
  section_name = 'NodeData'
  def do_parse(self,section):
    self.node_data = []
    self.data_name = section.lines[1]
    self.step_percents = float(section.lines[3])
    count = int(section.lines[7])
    self.node_data = map(lambda l:[float(i) for i in l.split()[1:]],
      section.lines[8:])


class GmshElementData(GmshSectionParser):
  section_name = 'ElementData'
  def do_parse(self,section):
    self.element_data = []
    self.data_name = section.lines[1]
    self.step_percents = float(section.lines[3])
    count = int(section.lines[7])
    self.element_data = map(lambda l:[float(i) for i in l.split()[1:]],
      section.lines[8:])


class Gmsh(object):
  """
  Class encapsulating an access to gmsh file structure
  creates an array 'sections' of GmshSection object instances
  by given filename
  """
  # a list of classes which represents a recognized sections
  section_mapping_objects = (GmshMeshFormat,GmshNodes,
    GmshElements,GmshNodeData,GmshElementData,)
  # a dictionary with the section name as a key
  # and corresponding class as a value
  section_mapping = dict()
  # construct this mapping
  for ob in section_mapping_objects:
    section_mapping[ob.section_name] = ob

  def __init__(self,filename = ""):
    super(Gmsh,self).__init__()
    self.sections = []
    if len(filename):
      self.parse_file(filename)
  
  def parse_file(self,filename):
    # read contents
    f = open(filename,'r')
    contents = f.read()
    f.close()
    # parse sections
    section_strings = self.parse_string(contents)
    self.sections =map(lambda section:GmshSection(section),section_strings) 
    self.create_section_objects()

  def create_section_objects(self):
    # loop over all sections
    for section in self.sections:
      # find appropriate class for a section
      class_for_section = self.section_mapping[section.name]
      if class_for_section: # ok, class found
        # create a class instance for the section
        object_for_section = class_for_section()
        # parse a section into the class instances
        object_for_section.parse_section(section)
        # check if such attribue is already exist
        if hasattr(self,section.name):
          # check if such attribute is of object type
          section_attr = getattr(self,section.name)
          if type(section_attr) is class_for_section:
            # convert to an array
            setattr(self,section.name,[section_attr])
            # append new section to this array
            section_attr = getattr(self,section.name)
            section_attr.append(object_for_section)
          elif type(section_attr) is list:
            section_attr.append(object_for_section)
            
        else: # no such attribute 
          setattr(self,section.name,object_for_section)

  def parse_string(self,str):
    # Match string, see comments for explanations
    # MATCH_OLD = '\$(?!End)(?:(?!\$End\w+).)*\$End\w+'
    MATCH = r"""\$(?!End) # Find string which starts with '$' but not 
                          # followed by 'End'
                (?:(?!\$End\w+).)* # Any content except $End{something}
                \$End\w+ # ... and closing tag $End{something}"""
    # re.DOTALL means '.' is every char including newline
    # re.VERBOSE means to allow comments inside regexp
    R = re.compile(MATCH,re.DOTALL|re.VERBOSE)
    return R.findall(str)

def parse_file(filename):
  return Gmsh(filename)


