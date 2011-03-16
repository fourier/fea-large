/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef __XML_LOADER_H__
#define __XML_LOADER_H__

#include "defines.h"
#include "fea_solver.h"

#ifdef USE_EXPAT

/* All known XML tags */
typedef enum xml_format_tags_enum {
  UNKNOWN_TAG,
  TASK,
  MODEL,
  MODEL_PARAMETERS,
  SOLUTION,
  ELEMENT_TYPE,
  LINE_SEARCH,
  ARC_LENGTH,
  INPUT_DATA,
  GEOMETRY,
  NODES,
  NODE,
  ELEMENTS,
  ELEMENT,
  BOUNDARY_CONDITIONS,
  PRESCRIBED_DISPLACEMENTS,
  PRESC_NODE
} xml_format_tags;

/* An input data structure used in parser */
typedef struct parse_data_tag {
  fea_task *task;
  fea_solution_params *fea_params;
  nodes_array *nodes;
  elements_array *elements;
  presc_boundary_array *presc_boundary;
  xml_format_tags parent_tag;
  char* current_text;
  int current_size;
} parse_data;


/*************************************************************/
/* Declarations of functions used                            */

/* Case-insensitive string comparsion procedure */
int istrcmp(char *s1,char *s2);

/* Convert particular string to the XML tag enum */
xml_format_tags tagname_to_enum(const XML_Char* name);

/*
 * Remove leading and trailing whitespaces from the string,
 * allocating null-terminated string as a result
 */
char *trim_whitespaces(const char* string,size_t size);

/* Expat start tag handler */
void expat_start_tag_handler(void *userData,
                             const XML_Char *name,
                             const XML_Char **atts);

/* Expat End tag handler */
void expat_end_tag_handler(void *userData,
                           const XML_Char *name);


/*    
 * Expat handler for the text between tags
 * Since this function could be called several times for the current tag,
 * it is necessary to store text somewhere. We use parse_data->current_text
 * pointer and parse_data->current_size for these purposes
 */
void expat_text_handler(void *userData,
                        const XML_Char *s,
                        int len);

/* test if an attribute name is what expected(attribute_name)
 * and increase pointer to the next attribute if yes*/
BOOL check_attribute(const char* attribute_name, const XML_Char ***atts);

/* model tag handler */
void process_model_type(parse_data* data, const XML_Char **atts);

/* model-parameters tag handler */
void process_model_params(parse_data* data, const XML_Char **atts);

/* solution tag handler */
void process_solution(parse_data* data, const XML_Char **atts);

/* element-type tag handler */
void process_element_type(parse_data* data, const XML_Char **atts);

/* line-search tag handler */
void process_line_search(parse_data* data, const XML_Char **atts);

/* arc-length tag handler */
void process_arc_length(parse_data* data, const XML_Char **atts);

/* nodes tag handler */
void process_nodes(parse_data* data, const XML_Char **atts);

/* node tag handler */
void process_node(parse_data* data, const XML_Char **atts);

/* elements tag handler */
void process_elements(parse_data* data, const XML_Char **atts);

/* take the node id/position from the element attributes
 * like 'node1' or 'node10'
 * returns -1 in case of wrong attribute name
 * but not skip it in this case!
 */
int node_position_from_attr(const XML_Char ***atts);

/* element tag handler */
void process_element(parse_data* data, const XML_Char **atts);

/* prescribed-displacements tag handler */
void process_prescribed_displacements(parse_data* data, const XML_Char **atts);

/* node tag handler */
void process_prescribed_node(parse_data* data, const XML_Char **atts);

/* main function for loading data using expat XML processor */
BOOL expat_data_load(char *filename,
                     fea_task **task,
                     fea_solution_params **fea_params,
                     nodes_array **nodes,
                     elements_array **elements,
                     presc_boundary_array **presc_boundary);


#endif /* USE_EXPAT */

#endif /* __XML_LOADER_H__ */
