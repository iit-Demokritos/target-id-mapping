package gr.demokritos.transformations;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class CreateTargetMappings {
	
	public static void main(String [] args) {
		create(args);
		enrich(args);
		enrichTTDs(args);
		enrichDisGenetIds(args);
		
	}
	
	public static void enrichDisGenetIds(String []args) {
		try {
        	String dtiFolder = args[0];
        	int n=0;
        	
			String targetmap1 = dtiFolder+"target-mappings_upd.tsv";
			String targetmap2 = dtiFolder+"target-mappings_latest.tsv";
			
			
			BufferedReader br = new BufferedReader(new FileReader(targetmap1));
			PrintWriter writerDrugMap= new PrintWriter(targetmap2, "UTF-8");

			String line= br.readLine();
			writerDrugMap.println(line+"\tdisgenet_gene_id");
			while ((line = br.readLine()) != null ){
				String[] values = line.split("\t");
				String unid = values[2];
				String disgenid = getDisGenetMapping(unid, dtiFolder); //DisGenet_uniprot_mapping.tsv
				String addition = disgenid;
				if (!line.endsWith("\t"))
					addition="\t"+addition;
				writerDrugMap.println(line+addition);
			}
			br.close();
			writerDrugMap.close();
		} catch(Exception e) {
			e.printStackTrace();
		}
		
	}
	
	static String getDisGenetMapping(String unip, String ddiPath) {
		
   		String id="null";
		String line;
		String drugMappings = ddiPath+"DisGenet_uniprot_mapping.tsv";
		try {
			BufferedReader br = new BufferedReader(new FileReader(drugMappings));
	   		
			while ((line = br.readLine()) != null ){
				if(line.contains(unip)) {
			        String[] values = line.split("\t");
			        id = values[1];
			        break;
				}    
			}
			br.close();
		}catch(Exception e) {
			e.printStackTrace();
		}
		
		return id;

	}
	
	public static void enrichTTDs(String []args) {

    	try {
        	String dtiFolder = args[0];
        	int n=0;
        	
        	HashMap<String,TargetEntry> targetmap = new HashMap<String,TargetEntry>();
        	Collection<String> uniProtNames;
        	Set<String> geneNames = new HashSet<String> ();
        
        	String targetMapping_uniprot_all=dtiFolder+"protein_exps/TTD/P2-01-TTD_uniprot_all.txt";
        	String targetMappingFile4=dtiFolder+"protein_exps/TTD/P1-01-TTD_target_download.txt";

          	TargetEntry targetEntry = new TargetEntry();
        	//old mappings file - backup all previous entries in a hashmap
        	String oldMappingFile=dtiFolder+"target-mappings.tsv";
       		FileReader mappingsFileReader1 = new FileReader(oldMappingFile);
       		BufferedReader mappingsFileLines1 = new BufferedReader(mappingsFileReader1);
       		String mLine=mappingsFileLines1.readLine();
       		while ((mLine=mappingsFileLines1.readLine())!=null) {
       			targetEntry = new TargetEntry();
      			String[] line = mLine.split("\t");
      			targetEntry.cui = line[0];
      			targetEntry.uniprot_name = line[1];
      			
      			targetEntry.uniprot_cid = line[2];
      			targetEntry.gene_name = line[3];
      			for (String g:targetEntry.gene_name.split(","))
      				geneNames.add(g);
      			targetEntry.ensemblegene_id=line[4];
      			targetEntry.ttd_id = line[5];
      			targetmap.put(targetEntry.uniprot_name, targetEntry);
    		}
    		mappingsFileLines1.close();
    		mappingsFileReader1.close();
    		
    		uniProtNames = targetmap.keySet();//.add(targetEntry.uniprot_name);
        	//Open TTD_uniprot_all file and look for uniprot names 
        	FileReader interactionsReader2 = new FileReader(targetMapping_uniprot_all);
        	BufferedReader fileLines2 = new BufferedReader(interactionsReader2);
        	String fileline;
        	//ignore headers lines
        	int i=0;
        	while((fileline = fileLines2.readLine() ) != null) {
        		if (++i<24)
        			continue;
        		if (fileline.contains("UNIPROID")) {
        			String[] line = fileline.split("\t");
        			String Uniprot_name= line[2];
        			if (uniProtNames.contains(Uniprot_name)) {
        				String ttd_id = line[0];
        				TargetEntry targetentry = targetmap.get(Uniprot_name);
        				if (targetentry.ttd_id.equals("null"))
        					targetentry.ttd_id=ttd_id;
        				else if (targetentry.ttd_id.contains(ttd_id))
        					;
        				else
        					targetentry.ttd_id+=","+ttd_id;
        				targetmap.put(Uniprot_name, targetentry);
        			}
             	}
        	}
           	fileLines2.close();
        	interactionsReader2.close();
           	
        	//Open TTD_target_download file and look for uniprot names or genes and synonyms
        	FileReader interactionsReader = new FileReader(targetMappingFile4);
        	BufferedReader fileLines = new BufferedReader(interactionsReader);
           	//ignore headers lines
        	i=0;
        	while((fileline = fileLines.readLine() ) != null) {
        		if (++i<41)
        			continue;
        		if (fileline.contains("UNIPROID")) {
         			String[] line = fileline.split("\t");
        			String Uniprot_name= line[2];
        			if (uniProtNames.contains(Uniprot_name)) {
        				String ttd_id = line[0];
        				TargetEntry targetentry = targetmap.get(Uniprot_name);
        				if (targetentry.ttd_id.equals("null"))
        					targetentry.ttd_id=ttd_id;
        				else if (targetentry.ttd_id.contains(ttd_id))
        					;
        				else
        					targetentry.ttd_id+=","+ttd_id;
        				targetmap.put(Uniprot_name, targetentry);
        			}
        		}
        		else if (fileline.contains("GENENAME")) {
        			String[] line = fileline.split("\t");
        			String geneName= line[2];
        			String ttdId = line[0];
        			if (geneNames.contains(geneName)) {
        				for (TargetEntry t:targetmap.values()) {
        					if (t.gene_name.equals(geneName)) {
        						targetEntry=t;
        						break;
        					}
        				}
        				if (targetEntry.ttd_id.equals("null"))
        					targetEntry.ttd_id=ttdId;
        				else if (targetEntry.ttd_id.contains(ttdId))
        					;
        				else
        					targetEntry.ttd_id+=","+ttdId;
        				targetmap.put(targetEntry.uniprot_name, targetEntry);
        			}
        			else {  //new TTD - Gene name record
        				targetEntry = new TargetEntry();
        				targetEntry.gene_name = geneName;
        				targetEntry.ttd_id=ttdId;
        				targetEntry.cui="null";
        				targetEntry.uniprot_cid="null";
        				targetEntry.ensemblegene_id="null";
        				targetEntry.uniprot_name="null";
       					targetmap.put(""+(n++), targetEntry);
        				
        			}
        		}
        		//else if (fileline.contains("SYNONYMS")) {
        			//;
        		//}
        		//else ignore line
        	}
        	fileLines.close();
        	interactionsReader.close();
          	
        	//Now update mappings file!
        	String updMappingFile=dtiFolder+"target-mappings_upd.tsv";
        	PrintWriter targetMappingsTSV= new PrintWriter(updMappingFile, "UTF-8");
        	targetMappingsTSV.println("CUI\tUniprot_Name\tUniprot_id\tGene_name\tEnsembleGeneId\tTTD_id");
        	for (String c:targetmap.keySet()) {
        		targetEntry = targetmap.get(c);
            	targetMappingsTSV.println(targetEntry.cui+"\t"+targetEntry.uniprot_name+"\t"+targetEntry.uniprot_cid+"\t"+targetEntry.gene_name+"\t"+targetEntry.ensemblegene_id+"\t"+targetEntry.ttd_id);
        	}
        	targetMappingsTSV.close();
    	}
    	catch (Exception e) {
    		e.printStackTrace();
    	}
	}
	
	
	
	
	
	
	
	public static void enrich(String [] args) {		
    	try {
        	String dtiFolder = args[0];
        	int n=0; 
        	
    	  	//These are the Metathesaurus + opentarget files with target mappings
        	String targetMappingFile1=dtiFolder+"OpenTargets.ORG/MRSAT-CLEAN.RRF";
        	String targetMappingFile2=dtiFolder+"OpenTargets.ORG/tractability_buckets-2021-01-12.tsv";
        	HashMap<String,TargetEntry> targetmap = new HashMap<String,TargetEntry>();
        	
        	Set<String> ensemblIds = new HashSet<String> ();
        	Set<String> cuis_uniprotIds= new HashSet<String> ();
        	
        	TargetEntry targetEntry = new TargetEntry();
        	//old mappings file - backup all previous entries in a hashmap
        	String oldMappingFile=dtiFolder+"target-mappings-mini.tsv";
       		FileReader mappingsFileReader1 = new FileReader(oldMappingFile);
       		BufferedReader mappingsFileLines1 = new BufferedReader(mappingsFileReader1);
       		String mLine=mappingsFileLines1.readLine();
       		while ((mLine=mappingsFileLines1.readLine())!=null) {
       			targetEntry = new TargetEntry();
      			String[] line = mLine.split("\t");
      			targetEntry.cui = line[0];
      			targetEntry.uniprot_name = line[1];
      	//		uniProtNames.add(targetEntry.uniprot_name);
      			targetEntry.uniprot_cid = line[2];
      			cuis_uniprotIds.add(targetEntry.uniprot_cid);
      			targetEntry.gene_name = line[3];
      			targetEntry.ensemblegene_id="null";
      			targetEntry.ttd_id = line[4];
      			targetmap.put(targetEntry.cui, targetEntry);
    		}
    		mappingsFileLines1.close();
    		mappingsFileReader1.close();
    		
    		//enrich with ensembleIds
      		FileReader mappingsFileReader2 = new FileReader(targetMappingFile1);
       		BufferedReader mappingsFileLines2 = new BufferedReader(mappingsFileReader2);
       		while ((mLine=mappingsFileLines2.readLine())!=null) {
       			if (mLine.contains("ENSEMBLGENE_ID")) {
       				String []line =  mLine.split("\\|");
       				String cui = line[0];
       				String ensemblegene_id = line[10];
       				targetEntry = targetmap.get(cui);
       				if (targetEntry==null) {
       					targetEntry = new TargetEntry();
       					targetEntry.uniprot_cid="null";
       				}	
       				//else {
       				//	targetEntry.ensemblegene_id = ensemblegene_id;
       				//	targetmap.put(targetEntry.uniprot_cid, targetEntry);
       				//}
       				targetEntry.cui=cui;
       				targetEntry.gene_name="null";
       				targetEntry.ensemblegene_id = ensemblegene_id;
       				ensemblIds.add(ensemblegene_id);
       				targetmap.put(cui, targetEntry);
       			}
    		}
    		mappingsFileLines2.close();
    		mappingsFileReader2.close();
       			
    		//enrich with more ensembleIds and gene names
      		FileReader mappingsFileReader3 = new FileReader(targetMappingFile2);
       		BufferedReader mappingsFileLines3 = new BufferedReader(mappingsFileReader3);
       		mLine=mappingsFileLines3.readLine();//header
       		while ((mLine=mappingsFileLines3.readLine())!=null) {
       			String []line = mLine.split("\t");
       			String ensemblegene_id = line[0];
       			String geneName = line[1];
       			String uniprotId = line[2];
   				geneName = geneName.replace("; ", ",");
   				if (ensemblIds.contains(ensemblegene_id)) {
       				for (TargetEntry t:targetmap.values()) {
       					if (t.ensemblegene_id.equals(ensemblegene_id)) {
       						targetEntry=t;
       						break;
       					}
       					//targetEntry = targetmapTMP.get(ensemblegene_id);
       				}	
       				targetEntry.gene_name=geneName;
       				if (!targetEntry.uniprot_cid.contains(uniprotId)) {
       					if (targetEntry.uniprot_cid.equals("null"))
       						targetEntry.uniprot_cid=uniprotId;
       					else 
       						targetEntry.uniprot_cid+=","+uniprotId;
       				}
       				targetmap.put(targetEntry.cui,targetEntry);
   				}
   				else if (cuis_uniprotIds.contains(uniprotId)) {
   					for (TargetEntry t:targetmap.values()) {
       					if (t.uniprot_cid.equals(uniprotId)) {
       						targetEntry=t;
       						break;
       					}
       				}
   					//targetEntry = targetmap.get(uniprotId);
					targetEntry.ensemblegene_id=ensemblegene_id;
					targetEntry.gene_name=geneName;
					targetmap.put(targetEntry.cui, targetEntry);
   				}		
   				else { //new record!
   					targetEntry = new TargetEntry();
   					targetEntry.ensemblegene_id=ensemblegene_id;
   					targetEntry.gene_name = geneName;
   					targetEntry.uniprot_cid=uniprotId;
   					targetEntry.cui="null";
   					targetmap.put(""+(n++), targetEntry);
   				}    		
   			}
    		mappingsFileLines3.close();
    		mappingsFileReader3.close();
 
        	
    		//Now update mappings file!
    		String updMappingFile=dtiFolder+"target-mappings.tsv";
    		PrintWriter targetMappingsTSV= new PrintWriter(updMappingFile, "UTF-8");
    		targetMappingsTSV.println("CUI\tUniprot_Name\tUniprot_id\tGene_name\tEnsembleGeneId\tTTD_id");
    		for (String c:targetmap.keySet()) {
    			targetEntry = targetmap.get(c);
        		targetMappingsTSV.println(targetEntry.cui+"\t"+targetEntry.uniprot_name+"\t"+targetEntry.uniprot_cid+"\t"+targetEntry.gene_name+"\t"+targetEntry.ensemblegene_id+"\t"+targetEntry.ttd_id);
    		}
    		targetMappingsTSV.close();
        	
    		
    	}catch (Exception e) {
    		e.printStackTrace();
    	}
	}
	
	public static void create(String [] args) {		
    	try {
    		 
        	String dtiFolder = args[0];
        	//This is the TTD file with drug-protein correlations
        	
       		HashMap<String,TargetEntry> proteinmap = new HashMap<String,TargetEntry>();
    		
    		//get mapping of all targets
    		HashMap<String,String> target_CUIs_Uniprot_ids = new HashMap<String,String>();
 
    	  	//This is the Metathesaurus file with target mappings
        	String targetMappingFile=dtiFolder+"OpenTargets.ORG/MRSAT-CLEAN.RRF";
        	String uniprotId, cui;
        	
    		FileReader mappingsFileReader = new FileReader(targetMappingFile);
    		BufferedReader mappingsFileLines = new BufferedReader(mappingsFileReader);
    		String mLine="";
    		int l=0;
    		while ((mLine=mappingsFileLines.readLine())!=null) {
    			l++;
   				if (mLine.contains("SWP")) {
    					String[] line = mLine.split("\\|");
    					cui = line[0];
    					if (mLine.contains("\"")) {   //many uniprotIds
    						line[10]=line[10].replaceAll("\"","");
    						line[10]=line[10].replaceAll("&#x7C;", ",");
    					}
    					uniprotId= line[10];
    					//System.out.println("SAVE "+ cui+" line[1]="+line[1]+" line[2]="+line[2]);
    					target_CUIs_Uniprot_ids.put(cui, uniprotId);
   				}	
    		}
    		mappingsFileLines.close();
    		mappingsFileReader.close();
    
        	System.out.println(l+" Rows examined. Mapped CUIs to Uniprot ids: "+target_CUIs_Uniprot_ids.size());
        	HashMap<String,String> uniprot_ids_names = getReverseUniProtmapping(dtiFolder, target_CUIs_Uniprot_ids.values());
        	System.out.println("Mapped uniprotIds to Uniprot names: "+uniprot_ids_names.size());
        	proteinmap = getTargetsMapObject(uniprot_ids_names, target_CUIs_Uniprot_ids, dtiFolder);
        	System.out.println("Mapped CUIs to TTD: "+proteinmap.size()); 
    		//all drug-target pairs - must divide to POS/NEG based on TTD groundtruth
   
        	//Create output TSV file
    		PrintWriter targetMappingsTSV= new PrintWriter(dtiFolder+"target-mappings-mini.tsv", "UTF-8");
    		targetMappingsTSV.println("CUI\tUniprot_Name\tUniprot_cid\tGene_name\tTTD_id");
    		for (String c:proteinmap.keySet()) {
        		TargetEntry targetEntry= proteinmap.get(c);
        		targetMappingsTSV.println(targetEntry.cui+"\t"+targetEntry.uniprot_name+"\t"+targetEntry.uniprot_cid+"\t"+targetEntry.gene_name+"\t"+targetEntry.ttd_id);
        	}
        	targetMappingsTSV.close();

    		
        }catch(Exception e) {
        	e.printStackTrace();
        }
	}
	

	//search for pair in ttdFile
	public static boolean searchPairInTTDFile(String ttdFile, String ttdIds) throws Exception{
		
		FileReader ttdFileReader = new FileReader(ttdFile);
		BufferedReader ttdFileLines = new BufferedReader(ttdFileReader);
		String line=ttdFileLines.readLine();
		while ((line=ttdFileLines.readLine())!=null) {
				if (line.contains(ttdIds)) {
					ttdFileLines.close();
					ttdFileReader.close();
					return true;
				}	
		}
		//System.out.println("Error! Cannot find TTD mapping for CUI: "+drugCUI);
		ttdFileLines.close();
		ttdFileReader.close();
		return false;
	}	
  		
	
	//search for the ttdId of a drug CUI
	public static String getMapping(String drugCUI, String file) throws Exception{

		FileReader mappingsFileReader = new FileReader(file);
		BufferedReader mappingsFileLines = new BufferedReader(mappingsFileReader);
		String mLine=mappingsFileLines.readLine();
		while ((mLine=mappingsFileLines.readLine())!=null) {
				if (mLine.contains(drugCUI)) {
					String []line  = mLine.split("\t");
					mappingsFileLines.close();
					mappingsFileReader.close();
					//System.out.println("Found in mappings file, CUI: "+drugCUI+", return TTD-id="+line[2]);
					return line[2];
				}	
		}
		//System.out.println("Error! Cannot find TTD mapping for CUI: "+drugCUI);
		mappingsFileLines.close();
		mappingsFileReader.close();
		return null;
	}



	public static HashMap<String, TargetEntry> getTargetsMapObject(HashMap<String,String> uniprot_ids_names, HashMap<String,String> target_CUIs_Uniprot_ids, String dtiFolder) throws Exception{
	
		//Save final mapping of ids in a Hashmap
		HashMap<String,TargetEntry> proteinmap = new HashMap<String,TargetEntry>();
		//This is the TTD file with target info
    	String targetInfoFile=dtiFolder+"protein_exps/TTD/TTD_target_info.txt";
        //another mapping
    	HashMap<String,String> uniprotNames_ttdIds = new HashMap<String,String>();
    	
		Collection<String> uniProtNames = uniprot_ids_names.values();
		
		System.out.println("Look for ttd ids for: "+uniProtNames.size()+" names");
		
		//Open TTD file and look for uniprot names 
		FileReader interactionsReader = new FileReader(targetInfoFile);
		BufferedReader fileLines = new BufferedReader(interactionsReader);
       	String fileline;//ignore headers line
    	while((fileline = fileLines.readLine() ) != null) {
    		if (fileline.contains("UNIPROID")) {
    	        String[] line = fileline.split("\t");
    	        if (line.length<3)
    	        	continue;
    	        String Uniprot_name= line[2];
    	        if (uniProtNames.contains(Uniprot_name)) {
    	        	String ttd_id = line[0];
    	        	uniprotNames_ttdIds.put(Uniprot_name, ttd_id);
    	        }	
    	    }
    	}
    	fileLines.close();
    	interactionsReader.close();
    	
    	System.out.println("Found ttd ids: "+uniprotNames_ttdIds.size());
    	System.out.println("Now adding mappings to an object");
		for (String cui: target_CUIs_Uniprot_ids.keySet()) {
			String uniprotId =  target_CUIs_Uniprot_ids.get(cui);
			String uniprotName=null;
			if (uniprotId.contains(",")) {
				String [] ids = uniprotId.split(",");
				for (String distinct_id:ids)
					if ((uniprotName = uniprot_ids_names.get(distinct_id))!=null)
						break;
			}	
			else
				uniprotName = uniprot_ids_names.get(uniprotId);
			String ttdId = null;
			if (uniprotName!=null) 
				ttdId = uniprotNames_ttdIds.get(uniprotName);
			TargetEntry t = new TargetEntry(); 
			t.cui=cui;
			t.gene_name=null;
			t.uniprot_name=uniprotName;
			t.ttd_id=ttdId;
			t.uniprot_cid=uniprotId;
			proteinmap.put(cui, t);
		}
		return proteinmap;
	}

	public static HashMap<String,String> getReverseUniProtmapping (String dtiFolder, Collection<String> uniProtIds) throws Exception{
	    
		String uniprotFile = dtiFolder+"protein_exps/TTD/uniprotmap";
		HashMap<String,String> uniprot_ids_names = new HashMap<String,String>();
		
		//split multiple uniprot ids
		Collection<String> upIds = new ArrayList<String>();
		upIds.addAll(uniProtIds);
		for (String id_entry:uniProtIds) 
			if(id_entry.contains(",")) {
				//System.out.println("multiple ids: "+id_entry);
				String [] ids = id_entry.split(",");
				for (String id:ids)
					upIds.add(id);
				upIds.remove(id_entry);
			}
		uniProtIds = new ArrayList<String>();
		uniProtIds.addAll(upIds);
		
		//Open mappings file
		FileReader uniprotReader = new FileReader(uniprotFile);
		BufferedReader inputLines = new BufferedReader(uniprotReader);
       	String inputLine; 
       	//int n=0; int g=0;
    	while((inputLine = inputLines.readLine() ) != null) {
    		//n++;
	       	//if(n==50000)
	       		//break;
	       	String[] line = inputLine.split("\t");  
    	    if(uniProtIds.contains(line[0])) {
    	    	//g++;
    	    	//if (g%1000==0)
    	    	//	System.out.println("found 1000 more id-name pairs, saving:"+line[0]+","+line[1]); 
    	     	uniprot_ids_names.put(line[0], line[1]);
    	    }
    	}
    	inputLines.close();
    	uniprotReader.close();
    	
		return uniprot_ids_names;
	}	

}
