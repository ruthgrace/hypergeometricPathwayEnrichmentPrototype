/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pathwayprototype;

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
import java.util.ArrayList;
import java.util.Set;
import java.util.HashMap;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Collections;
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
    final double PVALUE_CUTOFF = 0.0000001; // you cou  ld turn this into something related to the number of pathways
    /**
     * @param args the command line arguments
     */
    public Pathwayprototype() {
        
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
            geneset.retainAll(allgenes);
            
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
        int samplesize = geneset.size(), markeditems, populationsize = allgenes.size();
        Hypergeometric h;
        testedpathways = new ArrayList<TestedPathway>();
        double p;
        while (pathwaynames.hasNext()) {
            pathwayname = (String) pathwaynames.next();
            //get pathway
            commonGenes = new HashSet<String>((ArrayList<String>) genesets.get(pathwayname));
            markeditems = commonGenes.size();
            commonGenes.retainAll(geneset);
            //see if there are common genes
            if (commonGenes.size() > 0) {
                try {
                    h = new Hypergeometric(samplesize, populationsize, markeditems);
                    p = h.cdf((double) commonGenes.size());
                    //store p-value from hypergeometric test
                    //if p is less than some threshhold
                    //color all genes involved in the pathway
                    testedpathways.add(new TestedPathway(p, pathwayname));
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
                    p = h.cdf((double) commonGenes.size());
                    this.markPathwayGenesInGPML(WIKIPATHWAYSFOLDER, commonGenes, pathwayname, PATHWAYOUTPUTFOLDER, PVALUE_CUTOFF);
                    testedpathways.add(new TestedPathway(p, pathwayname));
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
            writer.println("Gene symbol" + "\t" + "P-value");
            TestedPathway t;
            int size = testedpathways.size();
            for (int i = 0; i < size; i++) {
                t = testedpathways.get(i);
                writer.println(t.getName() + "\t" + t.getP());
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
                    writer.print("|"+nodes.item(0).getNodeValue());
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
    /*public void wikipathways2GMT(String pathwayfolder, String OUTFILE) {
        Pattern p,p2; Matcher m,m2;
        String PATHWAYNAMEREGEX = "^<Pathway[^\"]*\"[^\"]*\"[ ]*Name[^\"]*\"([^\"]*)\".*$";
        String PATHWAYDESCRIPTIONREGEX = "^.*Comment Source[ ]*=[ ]*\"WikiPathways-description\">([^<]*).*$";
        String GENESYMBOLREGEX = "^.*DataNode[ ]*TextLabel[ ]*=[ ]*\"([^\"]*)\".*Type[ ]*=[ ]*\"GeneProduct\".*$";
        String DATANODE = "DataNode";
        BufferedReader br;
        String line;
        int counter = 0;
        final File folder = new File(pathwayfolder);
        try {
                PrintWriter writer = new PrintWriter(GENESETFOLDER+OUTFILE, "UTF-8");

                for (final File fileEntry : folder.listFiles()) {
                        br = new BufferedReader(new FileReader(pathwayfolder+fileEntry.getName()));
                        line = br.readLine();
                        if (line != null) {

                                line = line.trim();
                                p = Pattern.compile(PATHWAYNAMEREGEX);
                                m = p.matcher(line);
                                while (!m.find() && line != null) {
                                        line = br.readLine();
                                        m = p.matcher(line);
                                }

                                //counter ensures uniqueness
                                writer.print(fileEntry.getName() + "\t" + m.group(1)+"_"+counter);
                                counter++;
                                p = Pattern.compile(PATHWAYDESCRIPTIONREGEX);
                                m = p.matcher(line);
                                p2 = Pattern.compile(DATANODE);
                                m2 = p2.matcher(line);
                                while (!m.find() && line!=null) {
                                        if (m2.find()) {
                                                break;
                                        }
                                        line = br.readLine();
                                        m = p.matcher(line);
                                        m2 = p2.matcher(line);
                                }

                                if (!m2.find(0)){
                                        writer.print("\t"+m.group(1));
                                }
                                else {
                                        writer.print("\tNo description.");
                                }
                                p = Pattern.compile(GENESYMBOLREGEX);
                                while (line!=null) {
                                        m = p.matcher(line);

                                        if (m.find()) {
                                                writer.print("\t"+m.group(1));
                                        }
                                        line = br.readLine();
                                }
                                writer.println();
                        }

                }
                writer.close();
        }
        catch (Exception e) {
                System.out.println(e.getMessage());
        }
    }*/
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