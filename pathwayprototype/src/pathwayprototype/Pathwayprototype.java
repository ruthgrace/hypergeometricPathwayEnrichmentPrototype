/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pathwayprototype;

import java.nio.file.Path;
import java.net.*;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import javax.xml.parsers.*;
import org.w3c.dom.*;
import org.xml.sax.InputSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.xpath.*;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.File;
import java.util.*;
import java.util.regex.*;
import java.io.PrintWriter;
import jsc.distributions.Hypergeometric;
/**
 *
 * @author ruthgrace
 */
public class Pathwayprototype {
    private ArrayList<String> genes;
    private ArrayList<String> effects;
    private ArrayList<TestedPathway> testedpathways;
    private Set<String> geneset;
    private Set<String> allgenes;
    private Set<String> effectset;
    private HashMap genesets;
    private ArrayList<String> positionMap;
    final String GENE_REGEX="HGVS=([^:^(]*)";
    final String EFFECT_REGEX = "EFFECT=([^;]*)";
    final String TESTFILE1 = "src/resources/annotatedVCF/WGS-001-03.gatk.snp.indel.jv.vcf";
    final String TESTFILE2 = "src/resources/annotatedVCF/WGS-002-03.gatk.snp.indel.jv.vcf";
    final String TESTFILE3 = "src/resources/annotatedVCF/WGS-003-03.gatk.snp.indel.jv.vcf";
    final String GENESETFOLDER = "src/resources/genesets/";
    final String OUTPUTFILE = "output/enriched_pathways_and_Pvalues.txt";
    final String PATHWAYOUTPUTFOLDER = "output/gpml_files/";
    final String WIKIPATHWAYSGMTFILE = "wikipathways.gmt";
    final String GENES_NOT_IN_GENESETS_FILE = "output/genes_not_in_genesets.txt";
    final String OUTPUTDIR = "output/";
    final String PATHWAYSOUTPUTDIR = "output/gpml_files/";
    final String HTMLFOLDER = "html_files/";
    final double PVALUE_CUTOFF = 0.0000001; // you cou  ld turn this into something related to the number of pathways
    /**
     * @param args the command line arguments
     */
    public Pathwayprototype() {
        positionMap = new ArrayList<String>();
    }
    public static void main(String[] args) {
        String WIKIPATHWAYSFOLDER = "src/resources/wikipathways/";
        String VCFFILE = "src/resources/annotatedVCF/WGS-001-03.gatk.snp.indel.jv.vcf";
        String GENESETFILE = "src/resources/genesets/wikipathways.gmt";
        Pathwayprototype pathwayobject = new Pathwayprototype();
        
        //run with wikipathways
        pathwayobject.hypergeometricWithWikiPathways(VCFFILE,WIKIPATHWAYSFOLDER);
        
        //run with GeneSets
        //String GENESETFILE = "src/resources/genesets/Human_GO_AllPathways_no_GO_iea_symbol.gmt";
        //pathwayobject.hypergeometricWithGeneSets(GENESETFILE, VCFFILE);
        
    }
    public void hypergeometricWithGeneSets(String GENESETFILE, String VcfFile) {
        File outputdir = new File(OUTPUTDIR);
        if (!outputdir.exists()) {
            outputdir.mkdir();
        }
        File gpmloutputdir = new File(PATHWAYSOUTPUTDIR);
        if (!gpmloutputdir.exists()) {
            gpmloutputdir.mkdir();
        }
        this.readGeneSet(GENESETFILE);
        this.readVCF(VcfFile, GENE_REGEX, EFFECT_REGEX, GENES_NOT_IN_GENESETS_FILE);
        this.hypergeometricTest();
        this.outputEnrichedGeneList(OUTPUTFILE);
        System.out.println("Sample size: "+this.geneset.size()+", Population size: "+this.allgenes.size());
    }
    public void hypergeometricWithWikiPathways(String VcfFile, String pathwayfolder) {
        File outputdir = new File(OUTPUTDIR);
        if (!outputdir.exists()) {
            outputdir.mkdir();
        }
        
        File gpmloutputdir = new File(PATHWAYSOUTPUTDIR);
        if (!gpmloutputdir.exists()) {
            gpmloutputdir.mkdir();
        }
        wikipathways2GMT(pathwayfolder, WIKIPATHWAYSGMTFILE);
        this.readWikiPathwayGeneSet(WIKIPATHWAYSGMTFILE);
        this.readVCF(VcfFile, GENE_REGEX, EFFECT_REGEX, GENES_NOT_IN_GENESETS_FILE);
        this.hypergeometricWikiPathwaysTest(pathwayfolder);
        this.outputEnrichedGeneList(OUTPUTFILE);
        this.makeGPMLnodesblack(PATHWAYOUTPUTFOLDER);
        try {
            this.displayPathwayHyperlinks(this.OUTPUTDIR+this.HTMLFOLDER, pathwayfolder);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Sample size: "+this.geneset.size()+", Population size: "+this.allgenes.size());
    }
    public void readWikiPathwayGeneSet(String filename) {
        //store genesets as hashmap of hashset, hashmap is the pathway name.
        //store genes as hashmap of ArrayList filled with pathways.
        genesets = new HashMap<String,HashMap<String,Integer>>();
        allgenes = new HashSet<String>();
        String[] linecomponents = new String[0];
        try {
            BufferedReader br = new BufferedReader(new FileReader(this.GENESETFOLDER+filename));
            String line = br.readLine();
            String pathwayname;
            HashMap<String,Integer> pathwaygenes;
            while (line != null) {
                line = line.trim();
                //tab separated
                //caps name % humancyc % abbreviation, human readable name, gene symbols
                linecomponents = line.split("\t");
                if (linecomponents.length >2) {
                    
                    pathwayname = linecomponents[0] + "\t" + linecomponents[1];
                    pathwaygenes = new HashMap<String,Integer>();
                    int counter = 0;
                    for (int i = 2; i < linecomponents.length; i++) {
                        pathwaygenes.put(linecomponents[i],counter++);
                        allgenes.add(linecomponents[i]);
                    }
                    this.genesets.put(pathwayname,pathwaygenes);
                }
                line = br.readLine();
            }
        }
        catch (Exception e) {
            for (int counter = 0; counter < linecomponents.length; counter++) {
                System.out.println(linecomponents[counter]);
            }
            
            e.printStackTrace();
        }
    }
    public void readGeneSet(String filename) {
        //store genesets as hashmap of hashset, hashmap is the pathway name.
        //store genes as hashmap of ArrayList filled with pathways.
        genesets = new HashMap<String,ArrayList>();
        allgenes = new HashSet<String>();
        String[] linecomponents;
        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            String line = br.readLine();
            String pathwayname;
            ArrayList pathwaygenes;
            while (line != null) {
                line = line.trim();
                //tab separated
                //caps name % humancyc % abbreviation, human readable name, gene symbols
                linecomponents = line.split("\t");
                pathwayname = linecomponents[0] + "\t" + linecomponents[1];
                pathwaygenes = new ArrayList<String>();
                for (int i = 2; i < linecomponents.length; i++) {
                    pathwaygenes.add(linecomponents[i]);
                    allgenes.add(linecomponents[i]);
                }
                this.genesets.put(pathwayname,pathwaygenes);
                line = br.readLine();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    public void makeGPMLnodesblack(String folderpath) {
        try {
            final File folder = new File(folderpath);
            for (final File fileEntry : folder.listFiles()) {
                Path path = Paths.get(folderpath+fileEntry.getName());
                Charset charset = StandardCharsets.UTF_8;
                String content = new String(Files.readAllBytes(path), charset);
                //REPLACE COLOR REGEX WITH COLOR BLACK
                content = content.replaceAll("foo", "bar");
                Files.write(path, content.getBytes(charset));
                
                fileEntry.getName();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    public void displayPathwayHyperlinks(String htmlFolder, String PATHWAYFOLDER) throws URISyntaxException {
        //extract nodes (Label = gene symbol, ID = something unique, use position, color = black unless in tested pathways)
        try {
            File outputdir = new File(htmlFolder);
            if (!outputdir.exists()) {
                outputdir.mkdir();
            }
        
            Document doc;
            XPath xpath;
            PrintWriter writer;
            String xpathexpression, geneSymbol;
            NodeList nodes;
            Node node, graphics,point;
            String[] javascriptNodes=new String[1];
            String[] javascriptEdges=new String[1];
            
            NamedNodeMap attributes;
            final File folder = new File(PATHWAYFOLDER);
            String[] htmlFilePaths = new String[folder.listFiles().length];
            String[] pathwayNames = new String[folder.listFiles().length];
            String[] pathwayDescriptions = new String[folder.listFiles().length];
            int counter = 0;
            HashMap<String,String> nonGeneProductNodes;
            HashSet<String> geneProductNodes;
            
            for (final File fileEntry : folder.listFiles()) {
                System.out.println(fileEntry.getName());
                doc = DocumentBuilderFactory.newInstance().newDocumentBuilder().parse(new InputSource(PATHWAYFOLDER+fileEntry.getName()));
                xpath = XPathFactory.newInstance().newXPath();
                nonGeneProductNodes = new HashMap<String,String>();
                geneProductNodes = new HashSet<String>();
                xpathexpression = "/Pathway";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                
                if (nodes!=null && nodes.getLength() >0) {
                    attributes = nodes.item(0).getAttributes();
                    pathwayNames[counter]=attributes.getNamedItem("Name").getTextContent();
                }
                else {
                    System.out.println("could not find pathway name of "+fileEntry.getName());
                    pathwayNames[counter]="No Pathway Name";
                }
                htmlFilePaths[counter] = htmlFolder+fileEntry.getName().replaceFirst("\\.gpml","\\.html");
                
                xpathexpression = "/Pathway/Comment";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                //gene products
                if (nodes.getLength() >0 && nodes.item(0).getNodeType()==Node.ELEMENT_NODE){
                    int commentNodeCounter = 0;
                    Node commentNode = nodes.item(commentNodeCounter);
                    while (commentNode!=null && commentNode.getAttributes().getNamedItem("Source")!= null && !commentNode.getAttributes().getNamedItem("Source").getNodeValue().equals("WikiPathways-description")){
                        commentNode = nodes.item(++commentNodeCounter);
                    }
                    if (commentNode!=null) {
                        pathwayDescriptions[counter] = commentNode.getTextContent();
                    }
                    else {
                        pathwayDescriptions[counter] = "No Pathway Description";
                    }
                }
                
                counter++;
                
                xpathexpression = "/Pathway/DataNode";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                //gene products
                if (nodes.getLength() >0 && nodes.item(0).getNodeType()==Node.ELEMENT_NODE){
                    javascriptNodes = new String[nodes.getLength()];
                    for (int i = 0; i < javascriptNodes.length; i++) {
                        node = nodes.item(i);
                        javascriptNodes[i] = this.processNode(node, nonGeneProductNodes, geneProductNodes, this.positionMap);
                    }
                }
                
                //metabolites, etc (everything but gene products)
                xpathexpression = "/Pathway/Label";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                 if (nodes.getLength() >0 && nodes.item(0).getNodeType()==Node.ELEMENT_NODE){
                     for (int i = 0; i < javascriptNodes.length; i++) {
                        node = nodes.item(i);
                        if (node == null || !node.hasChildNodes()) {
                            continue;
                        }
                        this.processNode(node, nonGeneProductNodes, geneProductNodes, this.positionMap);
                     }
                     //remove trailing comma
                 }
                
                //interactions
                xpathexpression = "/Pathway/Interaction";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                //get nodes
                if (nodes.getLength() >0 && nodes.item(0).getNodeType()==Node.ELEMENT_NODE){
                    javascriptEdges = new String[nodes.getLength()];
                    for (int i = 0; i < javascriptEdges.length; i++) {
                        node = nodes.item(i);
                        graphics = node.getFirstChild();
                        while (graphics!=null && !graphics.getNodeName().equals("Graphics")) {
                            graphics = graphics.getNextSibling();
                        }
                        if (graphics == null ) {
                            System.out.println("edge did not have graphics component");
                            continue;
                        }
                        javascriptEdges[i] = this.processEdge(graphics, nonGeneProductNodes, geneProductNodes);
                    }
                }
                
                //shapes
                xpathexpression = "/Pathway/Shape";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                 if (nodes.getLength() >0 && nodes.item(0).getNodeType()==Node.ELEMENT_NODE){
                     for (int i = 0; i < javascriptNodes.length; i++) {
                        node = nodes.item(i);
                        if (node == null || !node.hasChildNodes()) {
                            continue;
                        }
                        this.processNode(node, nonGeneProductNodes, geneProductNodes, this.positionMap);
                     }
                 }
                
                writeJS(pathwayNames[counter-1], pathwayDescriptions[counter-1], htmlFolder, fileEntry.getName().replaceFirst("\\.gpml","\\.js"), javascriptNodes, javascriptEdges, nonGeneProductNodes);
                writeHTML(fileEntry.getName(), htmlFolder, pathwayNames[counter-1], pathwayDescriptions[counter-1]);
            }
            this.writeCSS(htmlFolder+"cytoscape_javascript_prototype.css");
            this.writeShowMore(htmlFolder);
            PathwayPickerDialog picker = new PathwayPickerDialog(pathwayNames, htmlFilePaths);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    private void writeJS(String pathwayName, String pathwayDescription, String folder, String fileName, String[] javascriptNodes, String[] javascriptEdges, HashMap<String,String> nonGeneProductNodes) {
        try {
            String jsfilebeginning = "$('#cy').cytoscape({\n  layout: {\n  name: 'preset',\n  positions: {";
            String jsfilestyle = "}\n},\nstyle: cytoscape.stylesheet()\n    .selector('node')\n      .css({\n        'content': 'data(name)',\n        'shape': 'rectangle',\n        'text-valign': 'center',\n        'background-color': 'white',\n        'background-opacity': 1,\n        'color': 'black',\n        'font-size': 10,\n        'text-outline-width': 0,\n        'border-width': 2,\n        'border-color': 'black',\n        'border-opacity': 1,\n        'text-outline-color': '#888',\n        'z-index': 2\n      })\n      .selector('$node > node')\n      .css({\n        'padding-top': '4px',\n        'padding-left': '4px',\n        'padding-bottom': '4px',\n        'padding-right': '4px',\n        'z-index': 0\n      })\n      .selector('node[isInhibitor]')\n      .css({\n        'color': 'red',\n        'border-color': 'red'\n      })\n      .selector('node[isMetabolite]')\n      .css({\n        'color': 'blue',\n        'border-color': 'blue'\n      })\n      .selector('node[placementNode]')\n      .css({\n        'width': 1,\n        'height': 1,\n        'background-opacity': 1,\n        'background-color': 'black',\n        'border-opacity': 0\n      })\n      .selector('node[width]')\n      .css({\n        'width': 'data(width)',\n        'height': 'data(height)',\n      })\n      .selector('node[shape]')\n      .css({\n        'background-opacity': 0,\n        'background-color': 'white',\n        'shape': 'data(shape)',\n        'z-index': 0\n      })\n      .selector('node > $node')\n      .css({\n        'z-index': 0\n      })\n      .selector('node[isLabel]')\n      .css({\n      	'background-opacity': 0,\n      	'border-opacity': 0,\n        'z-index': 1\n      })\n      .selector('node[color]')\n      .css({\n        'border-width': 4,\n        'border-color': 'green'\n      })\n      .selector('edge')\n      .css({\n        'target-arrow-shape': 'triangle',\n        'z-index': 3,\n        'width': 2,\n        'line-color': 'black',\n        'target-arrow-color': 'black'\n      })\n    .selector('edge[noArrowHead]')\n      .css({\n        'target-arrow-shape': 'none'\n      })\n    .selector('edge[teeArrowHead]')\n      .css({\n        'target-arrow-shape': 'tee',\n        'line-color': 'red',\n        'target-arrow-color': 'red'\n      })\n    .selector('edge[veeArrowHead]')\n      .css({\n        'target-arrow-shape': 'vee'\n      })\n    .selector('edge[circularArrowHead]')\n      .css({\n        'target-arrow-shape': 'circle',\n        'target-arrow-fill': 'hollow'\n      })\n    .selector(':selected')\n      .css({\n        'background-color': 'black',\n        'line-color': 'black',\n        'target-arrow-color': 'black',\n        'source-arrow-color': 'black'\n      })\n    .selector('.faded')\n      .css({\n        'background-opacity': 0.25,\n        'text-opacity': 0\n      })\n    .selector('node#supergrandparent')\n      .css({\n        'background-opacity': 0,\n        'border-opacity': 0,\n        'border-color': 'white',\n        'z-index': -1\n      })\n    .selector('node#superparent')\n      .css({\n        'background-opacity': 0,\n        'border-opacity': 0,\n        'border-color': 'white',\n        'z-index': -1\n      }),\nelements: {\n    nodes: [\n    { data: { id: 'supergrandparent' } },\n    { data: { id: 'superparent', parent: 'supergrandparent'} },";
            String jsfilemiddle = "    ],\n    edges: [";
            String jsfileend = "    ]\n  },\n  \n  ready: function(){\n    window.cy = this;\n    \n    // giddy up...\n    \n    cy.elements().unselectify();\n    \n    cy.on('tap', 'node', function(e){\n      var node = e.cyTarget; \n      var neighborhood = node.neighborhood().add(node);\n      \n      cy.elements().addClass('faded');\n      neighborhood.removeClass('faded');\n    });\n    \n    cy.on('tap', function(e){\n      if( e.cyTarget === cy ){\n        cy.elements().removeClass('faded');\n      }\n    });\n  }\n});";
            PrintWriter writer = new PrintWriter(folder+fileName, "UTF-8");
            writer.println(jsfilebeginning);
            //print node position map
            Iterator positions = positionMap.iterator();
            if (positions.hasNext()) {
                writer.println(positions.next());
            }
            while (positions.hasNext()) {
                writer.println(", "+positions.next());
            }
            writer.println(jsfilestyle);
            writer.print(javascriptNodes[0]);
            for (int i = 1; i < javascriptNodes.length; i++) {
                writer.print(",\n"+javascriptNodes[i]);
            }
            
            Iterator<String> it = nonGeneProductNodes.values().iterator();
            while (it.hasNext()) {
                writer.print(",\n"+it.next());
            }
            writer.println();
            writer.println(jsfilemiddle);
            
            //remove trailing comma
            if (javascriptEdges[javascriptEdges.length-1].substring(javascriptEdges[javascriptEdges.length-1].length() - 1).equals(",")) {
                javascriptEdges[javascriptEdges.length-1] = javascriptEdges[javascriptEdges.length-1].substring(0,javascriptEdges[javascriptEdges.length-1].length() - 1);
                
            }
            
            for (int i = 0; i < javascriptEdges.length; i++) {
                writer.println(javascriptEdges[i]);
            }
            writer.println(jsfileend);
            writer.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    private void writeShowMore(String folder) {
        try {
            PrintWriter writer = new PrintWriter(folder+"showmore.js");
            writer.print("$(document).ready(function() {\n    $(\".show-more a\").each(function() {\n    var $link = $(this);\n    var $content = $link.parent().prev(\"div.text-content\");\n\n    console.log($link);\n\n    var visibleHeight = $content[0].clientHeight;\n    var actualHide = $content[0].scrollHeight - 1;\n\n    console.log(actualHide);\n    console.log(visibleHeight);\n\n    if (actualHide > visibleHeight) {\n        $link.show();\n    } else {\n        $link.hide();\n    }\n});\n\n$(\".show-more a\").on(\"click\", function() {\n    var $link = $(this);\n    var $content = $link.parent().prev(\"div.text-content\");\n    var linkText = $link.text();\nconsole.log(\"click!\");\n    $content.toggleClass(\"short-text, full-text\");\n\n    $link.text(getShowLinkText(linkText));\n\n    return false;\n});\n\nfunction getShowLinkText(currentText) {\n    var newText = '';\n\n    if (currentText.toUpperCase() === \"SHOW MORE\") {\n        newText = \"Show less\";\n    } else {\n        newText = \"Show more\";\n    }\n\n    return newText;\n}\n});\n\n\n");
            writer.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    private void writeHTML(String gpmlFileName, String folder, String pathwayName, String pathwayDescription) {
        try {
            PrintWriter writer = new PrintWriter(folder+gpmlFileName.replaceFirst("\\.gpml","\\.html"));
            
            
            writer.print("<!DOCTYPE html>\n<html>\n<head>\n<meta name=\"description\" content=\"[Pathway Display]\" />\n<script src=\"http://ajax.googleapis.com/ajax/libs/jquery/1/jquery.min.js\"></script>\n<meta charset=utf-8 />\n<title>");
            writer.print(pathwayName);
            writer.print("</title>\n  <script src=\"http://cytoscape.github.io/cytoscape.js/api/cytoscape.js-latest/cytoscape.min.js\"></script>\n<script src=\"showmore.js\"></script>\n  <link rel=\"stylesheet\" href=\"./cytoscape_javascript_prototype.css\" />\n\n<div class=\"text-container\">\n<h1>");
            writer.print(pathwayName+"</h1><div class=\"text-content short-text\">");
            writer.println(pathwayDescription);
            writer.print("</div>\n    <div class=\"show-more\">\n        <a href=\"#\">Show more</a>\n    </div>\n    </div>\n\n</head>\n\n<body>\n\n\n  <div id=\"cy\">\n<script type=\"text/javascript\" src=\"./");
            writer.print(gpmlFileName.replaceFirst("\\.gpml","\\.js"));
            writer.println("\"></script>\n  </div>\n</body>\n</html>");
            writer.close();
            
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    private void writeCSS(String fileName) {
        try {
            PrintWriter writer = new PrintWriter(fileName);
            writer.println("body { \n  font: 14px helvetica neue, helvetica, arial, sans-serif;\n}\n\ndiv.text-container {\n\n    text-align: center;\n    margin: 0 auto;\n    width: 75%;    \n\n}\n\n.text-content{\n    line-height: 1em;\n\n    text-align: left;\n\n}\n\n.short-text {\n    overflow: hidden;\n    height: 2em;\n}\n\n.full-text{\n    height: auto;\n}\n\nh1 {\n    font-size: 24px;   \n}\n\n.show-more {\n    padding: 10px 0;\n    text-align: center;\n}\n\n#cy {\n  height: 100%;\n  width: 100%;\n  position: absolute;\n}");
            writer.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    private String processNode(Node node, HashMap<String,String> nonGeneProductNodes, HashSet<String> geneProducts, ArrayList<String> positionMap) {
        String nodeTextStart = "      { data: { id: '";
        String nodeTextName = ", name: '";
        String singleQuote = "'";
        String colon = ":";
        String space = " ";
        String nodeTextPlacementNode = ", placementNode: true";
        String nodeTextColor = ", color: true";
        String nodeTextWeight = ", weight: ";
        String nodeTextHeight = ", height: ";
        String nodeTextEnd = "} }";
        String startCurlyBrace = "{", endCurlyBrace = "}";
        String x = "x: ", y = ", y: ";
        String width = ", width: ", height = ", height: ";

        
        NamedNodeMap attributes = node.getAttributes();
        String nodeText = "";
        String graphID = "", geneSymbol = "";
        String shapeAttribute = ", shape: ";
        String labelAttribute = ", isLabel: true";
        String inhibitorAttribute = ", isInhibitor: true";
        String metaboliteAttribute = ", isMetabolite: true";
        String parentAttribute = ", parent: '";
        boolean addToNonGeneProductNodes = false;
        
        
        //if shape is double lined, call processNode twice
        //once with slightly smaller width/height, once with slightly larger
        if (node.getNodeName().equals("Shape")) {
            if (node.hasChildNodes()) {
                String nullString = "null";
                Node childNode = node.getFirstChild();
                while (childNode!=null && !childNode.getNodeName().equals("Attribute")) {
                    childNode = childNode.getNextSibling();
                }
                if (childNode!=null && childNode.getAttributes().getNamedItem("Key")!=null && childNode.getAttributes().getNamedItem("Key").getNodeValue().equals("org.pathvisio.DoubleLineProperty")) {
                    Node largerCopy = node.cloneNode(true), smallerCopy = node.cloneNode(true);
                    //make unique graphIDs
                    if (node.getAttributes().getNamedItem("GraphId")!=null) {
                        //give smaller node a unique graphID
                        smallerCopy.getAttributes().getNamedItem("GraphId").setNodeValue(smallerCopy.getAttributes().getNamedItem("GraphId").getNodeValue()+"_smaller_outline");
                        //designated larger copy as the smaller copy's parent so that they move together
                        ((Element) smallerCopy).setAttribute("parentID",node.getAttributes().getNamedItem("GraphId").getNodeValue());
                        ((Element) largerCopy).setAttribute("parentID",nullString);
                        
                    }
                    
                    //remove attribute child so that there are not infinite nested double lines
                    this.removeChildAttribute(largerCopy);
                    this.removeChildAttribute(smallerCopy);
                    
                    Node largerGraphics = this.getShapeGraphicsNode(largerCopy);
                    Node smallerGraphics = this.getShapeGraphicsNode(smallerCopy);
                    
                    //compound nodes cannot have custom sizes - remove size attributes of largeCopy
                    ((Element) largerGraphics).removeAttribute("Height");
                    ((Element) largerGraphics).removeAttribute("Width");
                    //change width/height
                    ((Element) smallerGraphics).setAttribute("Height", (Double.parseDouble(((Element) smallerGraphics).getAttribute("Height"))-4.0)+"");
                    ((Element) smallerGraphics).setAttribute("Width", (Double.parseDouble(((Element) smallerGraphics).getAttribute("Width"))-4.0)+"");
                    
                    return this.processNode(largerCopy, nonGeneProductNodes, geneProducts, positionMap) + "\n" + this.processNode(smallerCopy, nonGeneProductNodes, geneProducts, positionMap);
                }
            }
        }
        
        if (attributes!=null && attributes.getLength() > 0 && attributes.getNamedItem("GraphId")!=null) {
            
            nodeText+=nodeTextStart;
            graphID = getGraphRefID(node, nonGeneProductNodes, geneProducts, positionMap).trim();
            if (attributes.getNamedItem("Type")==null || !attributes.getNamedItem("Type").getNodeValue().equals("GeneProduct")) {
                addToNonGeneProductNodes = true; 
            }
            else {
                geneProducts.add(graphID);
            }
            
            nodeText+=graphID+singleQuote;
            String nullString = "null";
            if (attributes.getNamedItem("parentID")==null) {
                String superparentString = "superparent'";
                nodeText+=parentAttribute+superparentString;
            }
            else if (!attributes.getNamedItem("parentID").getNodeValue().equals(nullString)) {
                nodeText+=parentAttribute+attributes.getNamedItem("parentID").getNodeValue()+singleQuote;
            }
            
            nodeText+=nodeTextName;
            
            if (attributes.getNamedItem("TextLabel")!=null) {
                geneSymbol = attributes.getNamedItem("TextLabel").getNodeValue().trim().replace("\n"," ");
                nodeText+=geneSymbol;
            }
            nodeText+= singleQuote;
            if (this.geneset.contains(geneSymbol)) {
                nodeText+=nodeTextColor;
            }
            
            
            if (node.getAttributes().getNamedItem("Type")!=null) {
                if (node.getAttributes().getNamedItem("TextLabel")!=null && node.getAttributes().getNamedItem("TextLabel").getNodeValue().toLowerCase().contains("inhibitor")) {
                    nodeText+=inhibitorAttribute;
                }
                else if (node.getAttributes().getNamedItem("Type").getNodeValue().equals("Metabolite")) {
                    nodeText+=metaboliteAttribute;
                }
                
            }
            
            Node graphics = node.getFirstChild();
            while (graphics!=null && !graphics.getNodeName().equals("Graphics")) {
                graphics = graphics.getNextSibling();
            }
            
            if (graphics == null ) {
                nodeText+=nodeTextEnd;
                return nodeText;
            }
            attributes = graphics.getAttributes();
            
            
            if (node.getNodeName().equals("Label")) {
                nodeText+=labelAttribute;
            }
            else if (node.getNodeName().equals("Shape")) {
                nodeText+=shapeAttribute;
                if (attributes.getNamedItem("ShapeType")!=null) {
                    String shapetype = attributes.getNamedItem("ShapeType").getNodeValue();
                    if (shapetype.equals("RoundedRectangle")) {
                        nodeText+="'roundrectangle'";
                    }
                    else if (shapetype.equals("Oval")) {
                        nodeText+="'ellipse'";
                    }
                    else {
                        nodeText+="'"+shapetype.toLowerCase().trim()+"'";
                    }
                }
                else {
                    nodeText+="'rectangle'";
                }
            }
            
            //width/height
            if (attributes.getNamedItem("Width")!=null && attributes.getNamedItem("Height")!=null) {
                String w = attributes.getNamedItem("Width").getNodeValue();
                String h = attributes.getNamedItem("Height").getNodeValue();
                nodeText+=width + w + height + h;
            }
            
            
            
            nodeText+=nodeTextEnd;
            
            //position
            if (attributes.getNamedItem("CenterX")!=null && attributes.getNamedItem("CenterY")!=null) {
                String xpos = attributes.getNamedItem("CenterX").getNodeValue();
                String ypos = attributes.getNamedItem("CenterY").getNodeValue();
                this.positionMap.add(singleQuote+graphID+singleQuote+colon+space+startCurlyBrace+x+xpos+y+ypos+endCurlyBrace);
            }
            else if (attributes.getNamedItem("X")!=null && attributes.getNamedItem("Y")!=null) {
                String xpos = attributes.getNamedItem("X").getNodeValue();
                String ypos = attributes.getNamedItem("Y").getNodeValue();
                this.positionMap.add(singleQuote+graphID+singleQuote+colon+space+startCurlyBrace+x+xpos+y+ypos+endCurlyBrace);
            }
            
            
        }
        
        
        if (addToNonGeneProductNodes) {
            nonGeneProductNodes.put(graphID, nodeText);
        }
        
        return nodeText;
    }
    private Node getShapeGraphicsNode(Node node) {
        Node childNode = node.getFirstChild();
        while (childNode!=null && !childNode.getNodeName().equals("Graphics")) {
            childNode = childNode.getNextSibling();
        }
        return childNode;
    }
    private void removeChildAttribute(Node node) {
        Node childNode = node.getFirstChild();
        while (childNode!=null && !childNode.getNodeName().equals("Attribute")) {
            childNode = childNode.getNextSibling();
        }
        node.removeChild(childNode);
    }
    private String processEdge(Node edgeNode, HashMap<String,String> nonGeneProductNodes, HashSet<String> geneProducts) {
        boolean arrowHead = true;
        Node target, source, nodeChild = edgeNode.getFirstChild();
        
        //arbitrarily assign points as Source and Target of edge
        while (!nodeChild.getNodeName().equals("Point")) {
            nodeChild = nodeChild.getNextSibling();
        }
        target = nodeChild;
        nodeChild = nodeChild.getNextSibling();
        while (!nodeChild.getNodeName().equals("Point")) {
            nodeChild = nodeChild.getNextSibling();
        }
        source = nodeChild;
        
        //correct assignment of points as Source and Target of edge - targets have arrowheads
        if (source.getAttributes().getNamedItem("ArrowHead")!=null) {
            Node temp = target;
            target = source;
            source = temp;
        }
        
        String sourcetext = "      { data: { source: '";
        String targettext = "', target: '";
        String arrowheadtext = "', noArrowHead: true";
        String singlequote = "'", endlinetext = " } },", linebreak = "\n";
        String edgeText="";
        
        Node anchor = edgeNode.getFirstChild();
        while (anchor!=null && !anchor.getNodeName().equals("Anchor")) {
            anchor = anchor.getNextSibling();
        }
        if (anchor!=null) {
            edgeText+=makeEdgeText(edgeNode,source, anchor, nonGeneProductNodes, geneProducts);
            edgeText+=linebreak;
            edgeText+=makeEdgeText(edgeNode,anchor, target, nonGeneProductNodes, geneProducts);
        }
        else {
            //check if elbow edge (perpendicular)
            if (edgeNode.getAttributes().getNamedItem("ConnectorType")!=null && edgeNode.getAttributes().getNamedItem("ConnectorType").getNodeValue().equals("Elbow") && this.getElbowCoordinates(source, target)!=null) {
                //give elbow positioning
                String elbowGraphID = this.getElbowCoordinates(source, target);
                String[] elbowCoordinates = elbowGraphID.split(" ");
                //make a node for the elbow anchor
                Element elbow = (Element) source.cloneNode(false);
                elbow.setAttribute("X",elbowCoordinates[0]);
                elbow.setAttribute("Y",elbowCoordinates[1]);
                elbow.removeAttribute("GraphId");
                elbow.removeAttribute("GraphRef");
                elbow.setAttribute("Elbow","true");
                elbow.setAttribute("TextLabel","");
                edgeText+=makeEdgeText(edgeNode,source, elbow, nonGeneProductNodes, geneProducts);
                edgeText+=linebreak;
                edgeText+=makeEdgeText(edgeNode,elbow, target, nonGeneProductNodes, geneProducts);
                //make sure elbow arms are red if part of an inhibition arrow
            }
            else {
                edgeText+=makeEdgeText(edgeNode,source, target, nonGeneProductNodes, geneProducts);
            }
        }
        
        return edgeText;
    }
    private String getElbowCoordinates(Node source, Node target) {
        String coordinates = "";
        if (target.getAttributes().getNamedItem("X")!= null) {
            coordinates+=target.getAttributes().getNamedItem("X").getNodeValue();
        }
        else if (target.getAttributes().getNamedItem("CenterX")!=null) {
            coordinates+=target.getAttributes().getNamedItem("CenterX").getNodeValue();
        }
        else {
            return null;
        }
        coordinates+=" ";
        if (source.getAttributes().getNamedItem("Y")!= null) {
            coordinates+=source.getAttributes().getNamedItem("Y").getNodeValue();
        }
        else if (source.getAttributes().getNamedItem("CenterY")!=null) {
            coordinates+=source.getAttributes().getNamedItem("CenterY").getNodeValue();
        }   
        else {
            return null;
        }
        return coordinates;
    }
    private String makeEdgeText(Node graphics, Node source, Node target, HashMap<String,String> nonGeneProductNodes, HashSet<String> geneProducts) {
        String sourcetext = "      { data: { source: '";
        String targettext = ", target: '";
        String labelText = ", label: '";
        String noArrowHeadText = ", noArrowHead: true";
        String circularArrowHeadText = ", circularArrowHead: true";
        String veeArrowHeadText = ", veeArrowHead: true";
        String triangleArrowHeadText = ", triangleArrowHead: true";
        String teeArrowHeadText = ", teeArrowHead: true";
        String singlequote = "'", endlinetext = " } },";
        String inhibitionString = "mim-inhibition";
        String catalysisString = "mim-catalysis";
        String bindingString = "mim-binding";
        String haystackText = ", haystackStyle: true";
        
        String edgeText = sourcetext;
        edgeText+= this.getGraphRefID(source, nonGeneProductNodes, geneProducts, positionMap)+singlequote;
        edgeText+=targettext;
        edgeText+= this.getGraphRefID(target, nonGeneProductNodes, geneProducts, positionMap)+singlequote;
        
        //MIM arrowhead conventions can be found here: http://discover.nci.nih.gov/mim/mapDesc.html
        
        if (target.getAttributes().getNamedItem("ArrowHead")==null) { //no arrowhead
            edgeText+=noArrowHeadText;
        }
        else {
            edgeText+=labelText;
            
            String arrowHeadType = target.getAttributes().getNamedItem("ArrowHead").getNodeValue();
            String[] arrowHeadLabel = arrowHeadType.split("-");
            edgeText+=arrowHeadLabel[arrowHeadLabel.length-1]+singlequote;
            
            if (arrowHeadType.equals(catalysisString)) {
                edgeText+=circularArrowHeadText;
            }
            else if (arrowHeadType.equals(inhibitionString)) {
                edgeText+=teeArrowHeadText;
            }
            else if (arrowHeadType.equals(bindingString)) {
                edgeText+=veeArrowHeadText;
            }
            else {
                edgeText+=triangleArrowHeadText;
            }
        }
        if (graphics.getAttributes().getNamedItem("ConnectorType")!=null && graphics.getAttributes().getNamedItem("ConnectorType").getNodeValue().equals("Elbow")) {
            edgeText+=haystackText;
        }
        edgeText+=endlinetext;
        return edgeText;
    }
    private String getGraphRefID(Node node, HashMap<String,String> nonGeneProductNodes, HashSet<String> geneProducts, ArrayList<String> positionMap) {
        String graphID;
        boolean newNonGeneProduct = false;
        if (node.getAttributes().getNamedItem("GraphRef")!=null) {
            graphID = node.getAttributes().getNamedItem("GraphRef").getNodeValue().trim();
            
        }
        else if (node.getAttributes().getNamedItem("GraphId")!=null) {
            graphID = node.getAttributes().getNamedItem("GraphId").getTextContent().trim();
        }
        else if (node.getAttributes().getNamedItem("X")!=null){
            graphID = node.getAttributes().getNamedItem("X").getNodeValue() + " " + node.getAttributes().getNamedItem("Y").getNodeValue().trim();
        }
        else if (node.getAttributes().getNamedItem("CenterX")!=null) {
            graphID = node.getAttributes().getNamedItem("X").getNodeValue() + " " + node.getAttributes().getNamedItem("Y").getNodeValue().trim();
        }
        else {
            graphID = node.getAttributes().getNamedItem("Position").getNodeValue();
        }
        //check centerx and figure out some defautl... maybe position
        String nodeText = "      { data: { id: '"+ graphID + "', name: ''";
        if (!geneProducts.contains(graphID) && !nonGeneProductNodes.containsKey(graphID)) {
            //check if there is positional information
            if (node.getAttributes().getNamedItem("X")!=null){
                positionMap.add("'"+graphID+"': {x: "+node.getAttributes().getNamedItem("X").getNodeValue()+", y: "+node.getAttributes().getNamedItem("Y").getNodeValue()+"}");
            }
            else if (node.getAttributes().getNamedItem("CenterX")!=null) {
                positionMap.add("'"+graphID+"': {x: "+node.getAttributes().getNamedItem("CenterX").getNodeValue()+", y: "+node.getAttributes().getNamedItem("CenterY").getNodeValue()+"}");
            }
            newNonGeneProduct = true;
            nodeText+=", placementNode: true";
        }
        if (node.getAttributes().getNamedItem("parentID")==null){
            nodeText+=", parent: 'superparent'";
        } 
        if (node.getAttributes().getNamedItem("Elbow")!=null) {
            nodeText+=", isElbowEdge: true";
        }
        nodeText+=" } }";
        if (newNonGeneProduct) {
            nonGeneProductNodes.put(graphID, nodeText);
        }
        return graphID;
    }
    
    public void readVCF(String filename, String GENE_REGEX, String EFFECT_REGEX, String GENES_NOT_IN_GENESETS_FILE) {
        genes = new ArrayList<String>();
        effects = new ArrayList<String>();
        geneset = new HashSet<String>();
        effectset = new HashSet<String> ();
        Pattern p;
        Matcher m;
        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            String line = br.readLine();
            String temp;
            while (line != null) {
                //skip header/preamble lines starting with #
                line = line.trim();
                if (line.charAt(0)!='#') {
                    //8th column has extra information separated by ;
                    p=Pattern.compile(GENE_REGEX);
                    m = p.matcher(line);
                    if (m.find()) {
                        genes.add(m.group(1));
                        geneset.add(m.group(1));
                        p=Pattern.compile(EFFECT_REGEX);
                        m = p.matcher(line);
                        m.find();
                        effects.add(m.group(1));
                        effectset.add(m.group(1));
                    }
                }
                line = br.readLine();
            }
            HashSet<String> genesNotInGeneSets = new HashSet<String>(geneset);
            genesNotInGeneSets.removeAll(allgenes);
            PrintWriter writer = new PrintWriter(GENES_NOT_IN_GENESETS_FILE, "UTF-8");
            writer.println("Gene symbol");
            Iterator<String> nonGenesetGenes = genesNotInGeneSets.iterator();
            while (nonGenesetGenes.hasNext()) {
                writer.println(nonGenesetGenes.next());
            }
            writer.close();
            geneset.retainAll(allgenes); //geneset is the set of mutated genes
            
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void hypergeometricTest() {
        //for each pathway
        HashSet<String> commonGenes;
        String pathwayname;
        Iterator pathwaynames = genesets.keySet().iterator();
        //marked items is the number of genes in VCF file
        //samplesize is number of genes in VCF file AND pathway in question
        //populationsize is all genes in genesets
        int markeditems = geneset.size(), samplesize, populationsize = allgenes.size();
        Hypergeometric h;
        testedpathways = new ArrayList<TestedPathway>();
        double p;
        while (pathwaynames.hasNext()) {
            pathwayname = (String) pathwaynames.next();
            //get pathway
            commonGenes = new HashSet<String>((ArrayList<String>) genesets.get(pathwayname));
            samplesize = commonGenes.size();
            commonGenes.retainAll(geneset);
            //see if there are common genes
            if (commonGenes.size() > 0) {
                try {
                    h = new Hypergeometric(samplesize, populationsize, markeditems);
                    p = h.cdf((double) commonGenes.size());
                    //store p-value from hypergeometric test
                    //if p is less than some threshhold
                    //color all genes involved in the pathway
                    testedpathways.add(new TestedPathway(p, pathwayname,commonGenes));
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        
        //comparators.sort by cdf
        Collections.sort(testedpathways, new PathwayEnrichmentProbabilityComparator());
        //comparators sort by pathway name (HUMAN READABLE)
        //Collections.sort(testedpathways, new PathwayNameComparator());
    }
    public void hypergeometricWikiPathwaysTest(String WIKIPATHWAYSFOLDER) {
        //for each pathway
        HashSet<String> commonGenes;
        String pathwayname;
        Iterator pathwaynames = genesets.keySet().iterator();
        int samplesize = geneset.size(), markeditems, populationsize = allgenes.size();
        Hypergeometric h;
        testedpathways = new ArrayList<TestedPathway>();
        double p;
        while (pathwaynames.hasNext()) {
            pathwayname = (String) pathwaynames.next();
            //get pathway
            commonGenes = new HashSet<String>( ( (HashMap<String,Integer>) genesets.get(pathwayname) ).keySet());
            markeditems = commonGenes.size();
            commonGenes.retainAll(geneset);
            //see if there are common genes
            if (commonGenes.size() > 0) {
                try {
                    h = new Hypergeometric(samplesize, populationsize, markeditems);
                    p = 1.0-h.cdf((double) commonGenes.size());
                    this.markPathwayGenesInGPML(WIKIPATHWAYSFOLDER, commonGenes, pathwayname, PATHWAYOUTPUTFOLDER, PVALUE_CUTOFF);
                    testedpathways.add(new TestedPathway(p, pathwayname,commonGenes));
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        
        //comparators.sort by cdf
        Collections.sort(testedpathways, new PathwayEnrichmentProbabilityComparator());
        //comparators sort by pathway name (HUMAN READABLE)
        //Collections.sort(testedpathways, new PathwayNameComparator());
    }
    public void outputEnrichedGeneList(String filename) {
        System.out.println("Writing pathways and p-values to "+filename);
        try {
            PrintWriter writer = new PrintWriter(filename, "UTF-8");
            writer.println("Gene symbol\tP-value\tGenes");
            TestedPathway t;
            HashSet<String> pathwayGenes;
            Iterator it;
            int size = testedpathways.size();
            for (int i = 0; i < size; i++) {
                t = testedpathways.get(i);
                writer.print(t.getName() + "\t" + t.getP());
                pathwayGenes=t.getGenes();
                it = pathwayGenes.iterator();
                while (it.hasNext()) {
                    writer.print("\t"+it.next());
                }
                writer.println();
            }
            writer.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    public void wikipathways2GMT(String pathwayfolder, String OUTFILE) {
        
        try {
            Document doc;
            XPath xpath;
            String xpathexpression;
            NodeList nodes;
            String geneNodeName;
            String[] geneNames;
            final File folder = new File(pathwayfolder);
            PrintWriter writer = new PrintWriter(GENESETFOLDER+OUTFILE, "UTF-8");
            for (final File fileEntry : folder.listFiles()) {
                doc = DocumentBuilderFactory.newInstance().newDocumentBuilder().parse(new InputSource(pathwayfolder+fileEntry.getName()));
                xpath = XPathFactory.newInstance().newXPath();
                
                writer.print(fileEntry.getName());
                
                xpathexpression = "/Pathway/@Name";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                if (nodes.getLength() >0 && nodes.item(0).getNodeType()==Node.ATTRIBUTE_NODE){
                    writer.print("\t"+((Attr) nodes.item(0)).getValue());
                }
                else {
                    System.out.println("unable to find pathway name in "+fileEntry.getName());
                }
                
                xpathexpression = "/Pathway/Comment[@Source='WikiPathways-description']";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                if (nodes.getLength() >0){
                    writer.print("|"+nodes.item(0).getTextContent().replaceAll("\\n"," "));
                }
                
                xpathexpression = "/Pathway/DataNode[@Type='GeneProduct']";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                for (int counter = 0; counter < nodes.getLength(); counter++) {
                    if (nodes.item(counter).getNodeType()==Node.ELEMENT_NODE) {
                        geneNodeName = ((Element) nodes.item(counter)).getAttribute("TextLabel");
                        geneNames = geneNodeName.split("[\\/]|\\s");
                        for (int counter2 = 0; counter2<geneNames.length;counter2++) {
                            if (geneNames[counter2].trim().length() != 0) {
                                writer.print("\t"+geneNames[counter2]);
                            }
                        }
                    }
                }
                writer.println();
            }
            writer.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    public void markPathwayGenesInGPML (String WIKIPATHWAYSFOLDER, HashSet<String> commonGenes, String pathwayName, String GPMLOUTPUTFOLDER, double PVALUE_CUTOFF) {
        Collections.sort(testedpathways, new PathwayEnrichmentProbabilityComparator());
        int i = 0;
        String pathwayFileName = pathwayName.split("\t")[0];
        String gene,geneMatch;
        try {
            
            Document doc = DocumentBuilderFactory.newInstance().newDocumentBuilder().parse(new InputSource(WIKIPATHWAYSFOLDER + pathwayFileName));

            XPath xpath = XPathFactory.newInstance().newXPath();
            
            Iterator enrichedGenes = commonGenes.iterator();
            
            
            while (enrichedGenes.hasNext()) {
                
                gene = (String) enrichedGenes.next();
                gene = gene.trim();
                String xpathexpression = "/Pathway/DataNode[contains(@TextLabel,'"+gene+"') and @Type='GeneProduct']";
                NodeList nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                boolean found = false;
                for (int counter = 0; counter < nodes.getLength(); counter++) {
                    if (nodes.item(counter).getNodeType()==Node.ELEMENT_NODE) {
                        //check what textlabel is (regex)
                        geneMatch = ((Element) nodes.item(counter)).getAttribute("TextLabel");
                        //get graphic element
                        
                        Pattern p = Pattern.compile("(^|\\s+|[/ '\"]+)"+gene+"(\\s+|[/ '\"]+|$)");
                        Matcher m = p.matcher(geneMatch);
                        if (m.find()) {
                            NodeList graphicComponents = ((Element) nodes.item(counter)).getElementsByTagName("Graphics");
                            ((Element) graphicComponents.item(0)).setAttribute("LineThickness","5.0");
                            found = true;
                        }
                        
                    }
                }
                if (!found) {
                    System.out.println("could not find "+gene+" in "+pathwayFileName + ", nodelength: "+nodes.getLength());
                    //also not fidning things with slashes
                }
            }
            Transformer xformer = TransformerFactory.newInstance().newTransformer();
            xformer.transform(new DOMSource(doc), new StreamResult(new File(GPMLOUTPUTFOLDER + pathwayFileName)));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
class PathwayEnrichmentProbabilityComparator implements Comparator<TestedPathway> {
    @Override
    public int compare(TestedPathway a, TestedPathway b) {
        if (a.getP() == b.getP()) {
            return 0;
        }
        if (a.getP() > b.getP()) {
            return 1;
        }
        return -1;
    }
    
}
class PathwayNameComparator implements Comparator<TestedPathway> {
    private final String NAME_REGEX = "\t([^\t])";
    @Override
    public int compare(TestedPathway a, TestedPathway b) {
        String readableNameA,readableNameB;
        Pattern p = Pattern.compile(NAME_REGEX);
        Matcher m = p.matcher(a.getName());
        m.find();
        readableNameA = m.group(1);
        m = p.matcher(b.getName());
        m.find();
        readableNameB = m.group(1);
        return readableNameA.compareToIgnoreCase(readableNameB);
    }
}