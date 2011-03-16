/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */

#include "defines.h"
#include "xml_loader.h"

#ifdef USE_EXPAT
#include <expat.h>

/*************************************************************/
/* Definition of functions used                              */

int istrcmp(char *s1,char *s2)
{
  /* case insensitive comparison */
  int d;
  for (;;) {
#ifdef ASCII_CTYPE
    if (!isascii(*s1) || !isascii(*s2))
      d = *s1 - *s2;
    else
#endif
      d = (tolower((unsigned char) *s1) - tolower((unsigned char)*s2));
    if ( d != 0 || *s1 == '\0' || *s2 == '\0' )
      return d;
    ++s1;
    ++s2;
  }
  /*NOTREACHED*/
}

xml_format_tags tagname_to_enum(const XML_Char* name)
{
  if (!istrcmp((char*)name,"TASK")) return TASK;
  if (!istrcmp((char*)name,"MODEL")) return MODEL;
  if (!istrcmp((char*)name,"MODEL-PARAMETERS")) return MODEL_PARAMETERS;
  if (!istrcmp((char*)name,"SOLUTION")) return SOLUTION;
  if (!istrcmp((char*)name,"ELEMENT-TYPE")) return ELEMENT_TYPE;
  if (!istrcmp((char*)name,"LINE-SEARCH")) return LINE_SEARCH;
  if (!istrcmp((char*)name,"ARC-LENGTH")) return ARC_LENGTH;
  if (!istrcmp((char*)name,"INPUT-DATA")) return INPUT_DATA;
  if (!istrcmp((char*)name,"GEOMETRY")) return GEOMETRY;
  if (!istrcmp((char*)name,"NODES")) return NODES;
  if (!istrcmp((char*)name,"NODE")) return NODE;
  if (!istrcmp((char*)name,"ELEMENTS")) return ELEMENTS;
  if (!istrcmp((char*)name,"ELEMENT")) return ELEMENT;
  if (!istrcmp((char*)name,"BOUNDARY-CONDITIONS"))
    return BOUNDARY_CONDITIONS;
  if (!istrcmp((char*)name,"PRESCRIBED-DISPLACEMENTS"))
    return PRESCRIBED_DISPLACEMENTS;
  if (!istrcmp((char*)name,"PRESC-NODE")) return PRESC_NODE;
  return UNKNOWN_TAG;
}

char *trim_whitespaces(const char* string,size_t size)
{
  const char* end = string+size;
  char* result = (char*)0;
  int not_ws_start = 0;
  int not_ws_end = 0;
  const char* ptr = string;
  /* find starting non-whitespace character */
  while( isspace(*ptr++) && size-- ) not_ws_start++;
  if (size != 0 || not_ws_start == 0)
  {
    ptr--;
    /* find trailing non-whitespace character */
    while(isspace(*--end) && end != ptr) not_ws_end++;
    size = end-ptr+1;
    result = (char*)malloc(size+1);
    memcpy(result,ptr,size);
    result[size] = '\0';
  }
  return result;
}

/*
 * Functions called from expat_start/end_tag_handler
 * when the tag is known
 */
void process_begin_tag(parse_data* data, int tag,const XML_Char **atts);
void process_end_tag(parse_data* data, int tag);

void expat_start_tag_handler(void *userData,
                             const XML_Char *name,
                             const XML_Char **atts)
{
  parse_data* data = (parse_data*)userData;
  int tag = tagname_to_enum(name);
  if(tag != UNKNOWN_TAG)
    process_begin_tag(data,tag,atts);
}

void expat_end_tag_handler(void *userData,
                           const XML_Char *name)
{
  parse_data* data = (parse_data*)userData;
  int tag = tagname_to_enum(name);

  if (tag != UNKNOWN_TAG)
    process_end_tag(data,tag);
  /* clear tag text data at tag close */
  if (data->current_text)
    free(data->current_text);
  data->current_text = (char*)0;
  data->current_size = 0;
}

void expat_text_handler(void *userData,
                        const XML_Char *s,
                        int len)
{
  parse_data* data = (parse_data*)userData;
  char *ptr;
  if (len)
  {
    if (!data->current_text)    /* allocate memory for the text in tag */
    {
      data->current_text = (char*)malloc(len);
      ptr = data->current_text;
    }
    else                        /* reallocate/widen memory alread allocated */
    {
      data->current_text = (char*)realloc(data->current_text,
                                          data->current_size+len);
      ptr = data->current_text;
      ptr += data->current_size;  
    }
    /* append text to the end of allocated/reallocad buffer */
    /* and increase size variable */
    memcpy(ptr,s,len);
    data->current_size += len;
  }
}


BOOL expat_data_load(char *filename,
                     fea_task **task,
                     fea_solution_params **fea_params,
                     nodes_array **nodes,
                     elements_array **elements,
                     presc_boundary_array **presc_boundary)
{
  BOOL result = FALSE;
  FILE* xml_document_file;
  XML_Parser parser;
  size_t file_size = 0;
  size_t read_bytes = 0;
  parse_data parse;
  enum XML_Status status;
  char *file_contents = (char*)0;

  /* Try to open file */
  if (!(xml_document_file = fopen(filename,"rt")))
  {
    printf("Error, could not open file %s\n",filename);
    return FALSE;
  }
  /* Determine file size */
  if (fseek(xml_document_file,0,SEEK_END))
  {
    printf("Error reading file %s\n",filename);
    return FALSE;
  }
  file_size = ftell(xml_document_file);
  /* rewind to the begin of file */
  fseek(xml_document_file,0,SEEK_SET);

  /* Create parser */
  parser = XML_ParserCreate(NULL);
  /* Set handlers */
  XML_SetElementHandler(parser, &expat_start_tag_handler,
                        &expat_end_tag_handler);
  XML_SetCharacterDataHandler(parser,expat_text_handler);

  /* initialize data */
  parse.current_text = (char*)0;
  parse.current_size = 0;

  /* read whole file */
  file_contents = (char*)malloc(file_size);
  read_bytes = fread(file_contents,1,file_size,xml_document_file);
  if (errno)
  {
    free(file_contents);
    return FALSE;
  }
  fclose(xml_document_file);

  /* allocate parse data */
  parse.task = new_fea_task();
  parse.fea_params = new_fea_solution_params();
  parse.nodes = new_nodes_array();
  parse.elements = new_elements_array();
  parse.presc_boundary = new_presc_boundary_array();
  parse.current_size = 0;
  parse.current_text = (char*)0;
  /* set user data */
  XML_SetUserData(parser,&parse);

  /* call parser */
  status = XML_Parse(parser,file_contents,(int)read_bytes,1);
  free(file_contents);
  XML_ParserFree(parser);
  
  *task = parse.task;
  *fea_params = parse.fea_params;
  *nodes = parse.nodes;
  *elements = parse.elements;
  *presc_boundary = parse.presc_boundary;

  
  result = TRUE;
  return result;
}


BOOL check_attribute(const char* attribute_name, const XML_Char ***atts)
{
  BOOL result = FALSE;
  if(!istrcmp((char*)(**atts),(char*)attribute_name))
  {
    (*atts)++;
    result = TRUE;
  }
  return result;
}

void process_model_type(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  for (; *atts; atts++ )
  {
    if (check_attribute("name",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      if (!istrcmp(text,"A5")) 
      {
        data->task->model.model = MODEL_A5;
        data->task->model.parameters_count = 2;
      }
      else if (!istrcmp(text,"COMPRESSIBLE_NEOHOOKEAN"))
      {
        data->task->model.model = MODEL_COMPRESSIBLE_NEOHOOKEAN;
        data->task->model.parameters_count = 2;
      }
      else
      {
        
        printf("unknown model type %s\n",text);
      }
      if(text)
        free(text);
    }
  }
}

/* model-parameters tag handler */
void process_model_params(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  int count = 0;
  for (; *atts && count < data->task->model.parameters_count; atts++ )
  {
    atts++;
    text = trim_whitespaces(*atts,strlen(*atts));
    data->task->model.parameters[count] = atof(text);
    if(text)
      free(text);
    count++;
  }
}

/* solution tag handler */
void process_solution(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  for (; *atts; atts++ )
  {
    if (check_attribute("modified-newton",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      data->task->modified_newton =
        (!istrcmp(text,"yes") || !istrcmp(text,"true"))? TRUE: FALSE;
      if (text)
        free(text);
    }
    else if (check_attribute("task-type",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      if (!istrcmp(text,"CARTESIAN3D"))
        data->task->type = CARTESIAN3D;
      if (text)
        free(text);
    }
    else if (check_attribute("load-increments-count",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      data->task->load_increments_count = atoi(text);
      if (text)
        free(text);
    }
    else if (check_attribute("desired-tolerance",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      data->task->desired_tolerance = atof(text);
      if (text)
        free(text);
    }
  }
  data->parent_tag = SOLUTION;
}

/* element-type tag handler */
void process_element_type(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  for (; *atts; atts++ )
  {
    if (check_attribute("name",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      if (!istrcmp(text,"TETRAHEDRA10"))
        data->task->ele_type = TETRAHEDRA10;
      if (text)
        free(text);
    }
    else if (check_attribute("nodes-count",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      data->fea_params->nodes_per_element = atoi(text);
      if (text)
        free(text);
    }
    else if (check_attribute("nodes-count",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      data->fea_params->nodes_per_element = atoi(text);
      if (text)
        free(text);
    }
    else if (check_attribute("gauss-nodes-count",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      data->fea_params->gauss_nodes_count = atoi(text);
      if (text)
        free(text);
    }
  }
}

/* line-search tag handler */
void process_line_search(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  for (; *atts; atts++ )
  {
    if (check_attribute("max",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      data->task->linesearch_max = atoi(text);
      if (text)
        free(text);
    }
  }
}

/* arc-length tag handler */
void process_arc_length(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  for (; *atts; atts++ )
  {
    if (check_attribute("max",&atts))
    {
      text = trim_whitespaces(*atts,strlen(*atts));
      data->task->arclength_max = atoi(text);
      if(text)
        free(text);
    }
  }
}

/* nodes tag handler */
void process_nodes(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  int i = 0;
  /* set parameters only when 'nodes' tag is a child of the 'geometry' tag */
  if ( data->parent_tag == GEOMETRY )
  {
    for (; *atts; atts++ )
    {
      if (check_attribute("count",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        data->nodes->nodes_count = atoi(text);
        /* allocate storage for nodes */
        data->nodes->nodes =
          (real**)malloc(data->nodes->nodes_count*sizeof(real*));
        for (; i < data->nodes->nodes_count; ++ i)
          data->nodes->nodes[i] = (real*)malloc(MAX_DOF*sizeof(real));
        if (text)
          free(text);
      }
    }
    /* set parent tag to 'nodes' to recoginze an appropriate 'node' tag */
    data->parent_tag = NODES;
  }
}

/* node tag handler */
void process_node(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  real dofs[MAX_DOF];
  int id = -1;
  if ( data->parent_tag == NODES )
  {
    for (; *atts; atts++ )
    {
      if (check_attribute("id",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        id = atoi(text);
        if (text) free(text);
      }
      else if (check_attribute("x",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        dofs[0] = atof(text);
        if (text) free(text);
      }
      else if (check_attribute("y",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        dofs[1] = atof(text);
        if (text) free(text);
      }
      else if (check_attribute("z",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        dofs[2] = atof(text);
        if (text) free(text);
      }
    }
    if (id != -1)
      memcpy(data->nodes->nodes[id],dofs,sizeof(dofs));
  }
}

/* elements tag handler */
void process_elements(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  int i = 0;
  /* set parameters only when 'nodes' tag is a child of the 'geometry' tag */
  if ( data->parent_tag == GEOMETRY )
  {
    for (; *atts; atts++ )
    {
      if (check_attribute("count",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        data->elements->elements_count = atoi(text);
        /* allocate storage for elements */
        data->elements->elements = 
          (int**)malloc(data->elements->elements_count*sizeof(int*));
        for (; i < data->elements->elements_count; ++ i)
          data->elements->elements[i] =
            (int*)malloc(data->fea_params->nodes_per_element*sizeof(int));
        if (text) free(text);
      }
    }
    /* set parent tag to 'ELEMENTS' to recoginze an appropriate 'ELEMENT' tag */
    data->parent_tag = ELEMENTS;
  }
}

int node_position_from_attr(const XML_Char ***atts)
{
  int result = -1;
  char* pos = strstr(**atts,"node");
  if (pos == **atts)
  {
    result = atoi(pos+strlen("node"))-1;
    (*atts)++;
  }
  return result;
}

void process_element(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  int pos = -1;
  int element_size_bytes = data->fea_params->nodes_per_element*sizeof(int);
  int* element = (int*)malloc(element_size_bytes);
  int id = -1;
  if ( data->parent_tag == ELEMENTS )
  {
    for (; *atts; atts++ )
    {
      if (check_attribute("id",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        id = atoi(text);
        if (text) free(text);
      }
      if (-1 != (pos = node_position_from_attr(&atts)))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        element[pos] = atoi(text); 
        if (text) free(text);
      }
    }
    if (id != -1)
      memcpy(data->elements->elements[id],element,element_size_bytes);
  }
  free(element);
}


/* prescribed-displacements tag handler */
void process_prescribed_displacements(parse_data* data, const XML_Char **atts)
{
  char* text = (char*)0;
  int size;
  /* set parameters only when 'nodes' tag is a child of the 'geometry' tag */
  if ( data->parent_tag == BOUNDARY_CONDITIONS )
  {
    for (; *atts; atts++ )
    {
      if (check_attribute("count",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        data->presc_boundary->prescribed_nodes_count = atoi(text);
        size = data->presc_boundary->prescribed_nodes_count;
        size = size*sizeof(prescibed_boundary_node);
        /* allocate storage for prescribed nodes */
        data->presc_boundary->prescribed_nodes =
          (prescibed_boundary_node*)malloc(size);
        if (text) free(text);
      }
    }
    /* set parent tag to 'nodes' to recoginze an appropriate 'node' tag */
    data->parent_tag = PRESCRIBED_DISPLACEMENTS;
  }
}

/* node tag handler */
void process_prescribed_node(parse_data* data, const XML_Char **atts)
{
  /* <presc-node id="1" node-id="10" x="0" y="0" z="0" type="7"/> */
  char* text = (char*)0;
  prescibed_boundary_node node;
  int id = -1;
  if ( data->parent_tag == PRESCRIBED_DISPLACEMENTS )
  {
    for (; *atts; atts++ )
    {
      if (check_attribute("id",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        id = atoi(text);
        if (text) free(text);
      }
      else if (check_attribute("node-id",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        node.node_number = atoi(text);
        if (text) free(text);
      }
      else if (check_attribute("x",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        node.values[0] = atof(text);
        if (text) free(text);
      }
      else if (check_attribute("y",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        node.values[1] = atof(text);
        if (text) free(text);
      }
      else if (check_attribute("z",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        node.values[2] = atof(text);
        if (text) free(text);
      }
      else if (check_attribute("type",&atts))
      {
        text = trim_whitespaces(*atts,strlen(*atts));
        /* TODO: add proper conversion */
        node.type= (presc_boundary_type)atoi(text);
        if (text) free(text);
      }
    }
    if (id != -1)
      memcpy(&data->presc_boundary->prescribed_nodes[id],&node,sizeof(node));
  }
}

void process_begin_tag(parse_data* data, int tag,const XML_Char **atts)
{
  switch(tag)
  {
  case TASK:
    break;
  case MODEL:
    process_model_type(data,atts);
    break;
  case MODEL_PARAMETERS:
    process_model_params(data,atts);
    break;
  case SOLUTION:
    process_solution(data,atts);
    break;
  case ELEMENT_TYPE:
    process_element_type(data,atts);
    break;
  case LINE_SEARCH:
    process_line_search(data,atts);
    break;
  case ARC_LENGTH:
    process_arc_length(data,atts);
    break;
  case INPUT_DATA:
    data->parent_tag = INPUT_DATA;
    break;
  case GEOMETRY:
    data->parent_tag = GEOMETRY;
    break;
  case NODES:
    process_nodes(data,atts);
    break;
  case NODE:
    process_node(data,atts);
    break;
  case ELEMENTS:
    process_elements(data,atts);
    break;
  case ELEMENT:
    process_element(data,atts);
    break;
  case BOUNDARY_CONDITIONS:
    data->parent_tag = BOUNDARY_CONDITIONS;
    break;
  case PRESCRIBED_DISPLACEMENTS:
    process_prescribed_displacements(data,atts);
    data->parent_tag = PRESCRIBED_DISPLACEMENTS;
    break;
  case PRESC_NODE:
    process_prescribed_node(data,atts);
    break;
  default:
    break;
  };
}

void process_end_tag(parse_data* data, int tag)
{
  switch(tag)
  {
  case NODE:
    break;
  case ELEMENT:
    break;
  case PRESC_NODE:
    break;
  case MODEL:
  case SOLUTION:
  case INPUT_DATA:    
    data->parent_tag = TASK;
    break;
  case MODEL_PARAMETERS:
    data->parent_tag = MODEL;
    break;
  case ELEMENT_TYPE:
  case LINE_SEARCH:
  case ARC_LENGTH:
    data->parent_tag = SOLUTION;
    break;
  case GEOMETRY:
  case BOUNDARY_CONDITIONS:
    data->parent_tag = INPUT_DATA;
    break;
  case NODES:
  case ELEMENTS:
    data->parent_tag = GEOMETRY;
    break;
  case TASK:
  default:
    data->parent_tag = UNKNOWN_TAG;
    break;
  };
}


#endif /* USE_EXPAT */


