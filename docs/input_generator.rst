Input Generator
===============

.. raw:: html

   <form id="input-form" onsubmit ="generateFile(event)">

       <label for="filename">Input Filename:</label>
       <input 
           type="text" 
           id="filename" 
           onblur="processFilename('filename')" 
           placeholder="Enter filename">
       <br>
       <small id="filenameHelp" style="color: #555; font-size: 0.9em; margin-top: 5px; display: block;">
           Enter the desired file name for the SparcleQC input file that will be created 
       </small>
       <br>


       <label for="pdb_file">PDB File:</label>
       <input type="text" id="pdb_file" name="pdb_file" placeholder="Enter PDB file name or path" style="margin-bottom: 7px;" required>
       <small id="pdbHelp" style="color: #555; font-size: 0.9em; margin-top: 5px; display: block;">
           Enter the path to the complex PDB (amber) or protein PDB (charmm)<br>
           The path should be either absolute or relative to where the input file will be executed<br>
           E.x., "example.pdb", "user/documents/example.pdb"
       </small>

       <br>

       <label for="cutoff_radius">Cutoff Radius (Å):</label>
       <input type="number" id="cutoff_radius" name="cutoff_radius" placeholder="Enter a floating-point number" step="any" min ="0" required style = "width: 27ch;">
       <small id="cutoffHelp" style="color: #555; font-size: 0.9em; margin-top: 5px; display: block;">
           Enter the cutoff radius (float) that determines which atoms are included in the QM region 
       </small>
       <br>

       <label for="seed_ligand">Seed:</label>
       <small id="seedHelp" style="color: #555; font-size: 0.9em; margin-top: 5px; display: block;">
           Choose if the QM region should be defined by distance to any atom in the ligand or to a specific atom<br>
           If choosing a specific atom, specify the file that contains that atom with the corresponding atom serial number (atomid)     
       </small>
       <div style="display: flex; align-items: center; margin-left: 30px; margin-top: 10px;">
           <label for="seed_ligand" style="margin-right: 10px;">Ligand:</label>
           <input type="checkbox" id="seed_ligand" name="seed_ligand" style="margin-right: 20px;">
           <label for="seed_id" style="margin-right: 10px;">Seed ID:</label>
           <input type="number" id="seed_id" name="seed_id" placeholder="Enter a number for seed ID">
           <label for="seed_file" style="margin-left: 20px; margin-right: 10px;">Seed File:</label>
           <input type="text" id="seed_file" name="seed_file" placeholder="Enter a seed file name">
       </div>

       <br>
       <label for="charge_scheme">Charge Scheme:</label>
       <select id="charge_scheme" name="charge_scheme">
           <option value="Z1">Z1</option>
           <option value="Z2">Z2</option>
           <option value="Z3">Z3</option>
           <option value="DZ1">DZ1</option>
           <option value="DZ2">DZ2</option>
           <option value="DZ3">DZ3</option>
           <option value="BRC">BRC</option>
           <option value="BRC2">BRC2</option>
           <option value="BRCD">BRCD</option>
       </select>
       <small id="schemeHelp" style="color: #555; font-size: 0.9em; margin-top: 5px; display: block;">
           Choose a charge scheme that determines how boundary MM charges are redistributed <br>
           Specifics on each charge scheme can be found in the <a href="user_guide.html" style="color: #007bff; text-decoration: none;">user guide</a>
       </small>
       <br>


       <label for="forcefield">Forcefield:</label>
       <small id="ffHelp" style="color: #555; font-size: 0.9em; margin-top: 5px; display: block;">
           If CHARMM is selected, charges must be precomputed using CHARMM-GUI<br>
           If Amber is selected, SparcleQC will obtain point charges automatically 
       </small>
       <div id="forcefields" style="margin-left: 30px; margin-top: 7px;">
          <div>
              <input type="checkbox" id="forcefield_charmm" name="forcefield" value="charmm" onclick="toggleExclusiveCheckbox('forcefield_charmm')">
              <label for="forcefield_charmm">CHARMM</label>
          </div>
          <small id="charmmHelp" style="color: #555; font-size: 0.9em; margin-top: 5px; display: none;">
              Enter the path to the CHARMM parameter files<br>  
              The path should be either absolute or relative to where the input file will be executed
          </small>
          <div id="charmm-options">
              <div style="margin-left: 30px;">
                  <label for="charmm_rtf" style="margin-right: 10px;">Topology Path:</label>
                  <input type="text" id="charmm_rtf" name="charmm_rtf" placeholder="Enter CHARMM RTF">
              </div>
              <div style="margin-left: 30px; margin-top: 10px;">
                  <label for="charmm_prm" style="margin-right: 10px;">Parameter Path:</label>
                  <input type="text" id="charmm_prm" name="charmm_prm" placeholder="Enter CHARMM PRM">
              </div>
          </div>
          <div style="margin-left: 0px;">
              <input type="checkbox" id="forcefield_amber" name="forcefield" value="amber" onclick="toggleExclusiveCheckbox('forcefield_amber')">
              <label for="forcefield_amber">Amber</label>
          </div>
          <small id="amberHelp" style="color: #555; font-size: 0.9em; margin-top: 5px; display: none;">
              Enter the desired Amber forcefield (e.g. ff19SB)<br>
              If other forcefields are needed to obtain point charges for the system enter them here <br>
              If the checkbox below is checked, SparcleQC will cap the terminal residues with ACE and NME  
          </small>
          <div id="amber-options" style="margin-top: 10px;">
              <div style="margin-left: 30px;">
                  <label for="amber_ff" style="margin-right: 10px;">Forcefield:</label>
                  <input type="text" id="amber_ff" name="amber_ff" placeholder="Enter Amber FF">
              </div>
              <div style="margin-left: 30px; margin-top: 10px;">
                  <label for="other_amber_ffs" style="margin-right: 10px;">Other Forcefields (Optional):</label>
                  <input type="text" id="other_amber_ffs" name="other_amber_ffs" placeholder="Enter other Amber FFs">
              </div>
              <div style="display: flex; align-items: center; margin-left: 30px;">
                  <label for="cap" style="margin-right: 10px;">Cap Terminal Residues?</label>
                  <input type="checkbox" id="precapbox" name="precap">
              </div>
          </div>
       </div>
       <br>
       
       <label for="water-charges-header" class="section-header">Water Charges:</label>
       <small id="waterHelp" style="color: #555; font-size: 0.9em; margin-top: 5px; display: block;">
           Enter the desired water model (e.g. OPC)<br>
           If you wish to override these charges, check the box below and you will be able to manually add charges for either a 3 or 4 point water
       </small>
       <div style="margin-left: 30px; margin-top: 7px;">
           <div style="margin-bottom: 10px;">
               <label for="water_model">Water Model:</label>
               <input type="text" id="water_model" name="water_model" placeholder="Enter water model" required>
           </div>
           <div style="margin-bottom: 10px;">
               <input type="checkbox" id="add_water_charges" onclick="toggleWater()">
               <label for="add_water_charges">Add Your Own Water Charges (Optional)</label>
           </div>
           <div style="display: none; margin-bottom: 10px;" id="tfield">
               <input type="checkbox" id="three_point_water" onclick="toggleWaterCharges(this)">
               <label for="three_point_water_model">3-point Water Model</label>
           </div>
           <div style="display: none; margin-bottom: 10px;" id="ffield">
               <input type="checkbox" id="four_point_water" onclick="toggleWaterCharges(this)">
               <label for="four_point_water_model">4-point Water Model</label>
           </div>
           <div style="display: none; margin-bottom: 20px;" id="o_charge_field">
               <label for="o_charge">Oxygen Charge:</label>
               <input type="number" id="o_charge" name="o_charge" placeholder="Enter oxygen charge" step="any">
           </div>
           <div style="display: none; margin-bottom: 20px;" id="h_charge_field">
               <label for="h_charge">Hydrogen Charge:</label>
               <input type="number" id="h_charge" name="h_charge" placeholder="Enter hydrogen charge" step="any">
           </div>
           <div style="display: none; margin-bottom: 20px;" id="ep_charge_field">
               <label for="ep_charge">Extra Point Charge:</label>
               <input type="number" id="ep_charge" name="ep_charge" placeholder="Enter extra point charge" step="any">
           </div>
       </div>

       <label for="software">Software:</label>
       <small id="softwareHelp" style="color: #555; font-size: 0.9em; margin-top: 5px; display: block;">
           Choose the desired quantum chemistry software<br>
           After choosing, optional software specific options will be available 
       </small>
       <div id="software" style="margin-left: 30px;">
           <label style="display: block; margin-top: 7px;">
               <input type="checkbox" id="software_nwchem" name="software" onclick="toggleSoftware('nwchem')"> NWChem
           </label>
           <div id="nwchem-options" style="display: none; margin-left: 20px;">
               <label style="margin-bottom: 10px; margin-top: 7px;" for="nwchemoptions">NWChem Settings (Optional):</label><br>
               <label for="nwchem_scratch">Scratch Directory:</label>
               <input type="text" id="nwchem_scratch" placeholder="Default: None" style="margin-bottom: 10px;"><br>
               
               <label for="nwchem_perm">Permanent Directory:</label>
               <input type="text" id="nwchem_perm" placeholder="Default: None" style="margin-bottom: 10px;"><br>
               
               <label for="nwchem_scf">SCF Options Dictionary:</label>
               <input type="text" id="nwchem_scf" placeholder="Default: None" style="margin-bottom: 10px;"><br>
               
               <label for="nwchem_dft">DFT Options Dictionary:</label>
               <input type="text" id="nwchem_dft" placeholder="Default: {'xc':'b3lyp'}" style="margin-bottom: 10px;"><br>
               
               <label for="mem">Memory:</label>
               <input type="text" id="mem" placeholder="Default: 32 GB" style="margin-bottom: 10px;">
           </div>
       
           <label style="display: block;">
               <input type="checkbox" id="software_qchem" name="software" onclick="toggleSoftware('qchem')"> Q-Chem
           </label>
           <div id="qchem-options" style="display: none; margin-left: 20px;">
               <label style="margin-bottom: 10px; margin-top: 7px;" for="qchemoptions">Q-Chem Settings (Optional):</label><br>
               <label for="qchem_options">Additional Options Dictionary:</label>
               <input type="text" size = "25" id="qchem_options" placeholder="Default: {'JOBTYPE': 'xsapt or sp'}" style="margin-bottom: 10px;"><br>
               <label for="qchem_sapt">SAPT Options (Dictionary):</label>
               <input type="text" id="qchem_sapt" size = "35" placeholder="Default: {} or {‘algorithm’:’ri-mo’,’basis’:’dimer’}" style="margin-bottom: 10px;">
           </div>
       
           <label style="display: block;">
               <input type="checkbox" id="software_psi4" name="software" onclick="toggleSoftware('psi4')"> Psi4
           </label>
           <div id="psi4-options" style="display: none; margin-left: 20px;">
               <label style="margin-bottom: 10px; margin-top: 7px;" for="psi4options">Psi4 Settings (Optional):</label><br>
               <label><input type="checkbox" id="fisapt_partition" style="margin-bottom: 10px;"> FISAPT Partition</label><br>
               <label><input type="checkbox" id="do_fsapt" style="margin-bottom: 10px;"> Do FSAPT</label><br>
               <label for="psi4options">Additional Options Dictionary:</label>
               <input type="text" id="psi4options" placeholder="Default: {}" style="margin-bottom: 10px;"><br>
               <label for="num_threads">Num Threads:</label>
               <input type="number" id="num_threads" placeholder="Default: 1" style="margin-bottom: 10px;"><br>
               <label for="memory">Memory:</label>
               <input type="text" id="memory" placeholder="Default: 32 GB" style="margin-bottom: 10px;">
           </div>
       </div>
       <br>
   
 

       <label for="ligand_charge">Ligand Charge:</label>
       <input type="number" id="ligand_charge" name="ligand_charge" placeholder="Enter ligand charge" step="1" required>
       <small id="ligHelp" style="color: #555; font-size: 0.9em; margin-top: 5px; display: block;">
           Enter the charge of the ligand 
       </small>
       <br>

       <label for="level_of_theory">Level of Theory:</label>
       <small id="theoryHelp" style="color: #555; font-size: 0.9em; margin-top: 5px; display: block;">
           Enter desired method (e.g. hf) and basis set (e.g. cc-pvdz) for the QM computation  
       </small>
       <div style="display: flex; align-items: center; margin-left: 30px; margin-top: 7px;">
           <label for="method" style="margin-right: 10px; margin-top: 7px;">Method:</label>
           <input type="text" id="method" name="method" placeholder="Enter method" required>
           <label for="basis_set" style="margin-left: 20px; margin-right: 10px;">Basis Set:</label>
           <input type="text" id="basis_set" name="basis_set" placeholder="Enter basis set" required>
       </div>
       <br>


       <label for="other_features">Other Features:</label>
       <small id="featureHelp" style="color: #555; font-size: 0.9em; margin-top: 5px; display: block;">
           Additional optional features are outlined in the <a href="user_guide.html" style="color: #007bff; text-decoration: none;">user guide</a> 
       </small>
       <div style="margin-left: 30px; margin-top: 7px;">
           <label for="template_path">Template Path (Optional):</label>
           <input type="text" id="template_path" size = "35" name="template_path" placeholder="Enter/Path/To/Template/cx_autocap_fixed.pdb">
       </div>
       <div style="margin-left: 30px; margin-top: 7px;">
           <label for="cp_correction">Counterpoise Correct?</label>
           <input type="checkbox" id="cp_correction">
       </div>
       <br>


       <button type="submit">Download</button>
   </form>




   <script>

       function toggleWater() {
           var addWaterCharges = document.getElementById("add_water_charges").checked;

           if (addWaterCharges) {
               document.getElementById("tfield").style.display = "inline-block";
               document.getElementById("ffield").style.display = "inline-block";
           } else {
               document.getElementById("tfield").style.display = "none";
               document.getElementById("ffield").style.display = "none";
               // Hide all charge fields
               document.getElementById("o_charge_field").style.display = "none";
               document.getElementById("h_charge_field").style.display = "none";
               document.getElementById("ep_charge_field").style.display = "none";
               const o_charge = document.getElementById("o_charge");
               const h_charge = document.getElementById("ep_charge");
               const ep_charge = document.getElementById("h_charge");
               o_charge.value = '';
               h_charge.value = '';
               ep_charge.value = '';
           }
       }

      


       function toggleWaterCharges(checkBox) {
           const threePointCheckbox = document.getElementById("three_point_water");
           const fourPointCheckbox = document.getElementById("four_point_water");

           if (checkBox === threePointCheckbox) {
               if (checkBox.checked) {
                   document.getElementById("o_charge_field").style.display = "block";
                   document.getElementById("h_charge_field").style.display = "block";
                   document.getElementById("ep_charge_field").style.display = "none"; 
                   fourPointCheckbox.checked = false; // Uncheck 4-point checkbox
                   const ep_charge = document.getElementById("ep_charge");
                   ep_charge.value = '';
               } else {
                   document.getElementById("o_charge_field").style.display = "none";
                   document.getElementById("h_charge_field").style.display = "none";
                   const o_charge = document.getElementById("o_charge");
                   const h_charge = document.getElementById("ep_charge");
                   o_charge.value = '';
                   h_charge.value = '';
               }
           } else if (checkBox === fourPointCheckbox) {
               if (checkBox.checked) {
                   document.getElementById("o_charge_field").style.display = "block";
                   document.getElementById("h_charge_field").style.display = "block";
                   document.getElementById("ep_charge_field").style.display = "block"; 
                   threePointCheckbox.checked = false; // Uncheck 3-point checkbox
               } else {
                   document.getElementById("o_charge_field").style.display = "none";
                   document.getElementById("h_charge_field").style.display = "none";
                   document.getElementById("ep_charge_field").style.display = "none";
                   const o_charge = document.getElementById("o_charge");
                   const h_charge = document.getElementById("ep_charge");
                   const ep_charge = document.getElementById("h_charge");
                   o_charge.value = '';
                   h_charge.value = '';
                   ep_charge.value = '';
               }
           }
       }

 
       


       function processFilename(inputId) {
           const inputField = document.getElementById(inputId);
           let filename = inputField.value.trim();
       
           if (!filename) {
               // If the input is empty, set the default filename
               filename = "output.txt";
           } else if (!filename.includes('.')) {
               // If there's no extension, add .txt
               filename += ".txt";
           }
       
           // Update the input field with the processed filename
           inputField.value = filename;
       }

       function toggleExclusiveCheckbox(selectedCheckboxId) {
           const forcefieldIds = ['charmm', 'amber'];
       
           forcefieldIds.forEach(forcefield => {
               const checkbox = document.getElementById(`forcefield_${forcefield}`);
               const optionsDiv = document.getElementById(`${forcefield}-options`);
               const notes = document.getElementById(`${forcefield}Help`);
       
               if (`forcefield_${forcefield}` === selectedCheckboxId) {
                   if (checkbox.checked) {
                       optionsDiv.style.display = "block"; // Show the selected options
                       notes.style.display = 'block';
                   } else {
                       optionsDiv.style.display = "none"; // Hide if unchecked
                       notes.style.display = "none"; // Hide if unchecked
                       const inputs = optionsDiv.querySelectorAll('input, select, textarea'); // Get all input elements
                       inputs.forEach(input => {
                           if (input.type === 'checkbox') {
                               input.checked = false; // Uncheck checkboxes
                           } else if (input.type === 'text') {
                               input.value = ''; // Clear textboxes
                           }
                       });

                   }
               } else {
                   const otherCheckbox = document.getElementById(`forcefield_${forcefield}`);
                   const otherOptionsDiv = document.getElementById(`${forcefield}-options`);
                   const notes = document.getElementById(`${forcefield}Help`);
                   otherCheckbox.checked = false; // Uncheck other checkboxes
                   otherOptionsDiv.style.display = "none"; // Hide other options
                   notes.style.display = "none"; // Hide other options
                   const inputs = optionsDiv.querySelectorAll('input, select, textarea'); // Get all input elements
                   inputs.forEach(input => {
                       if (input.type === 'checkbox') {
                           input.checked = false; // Uncheck checkboxes
                       } else if (input.type === 'text') {
                           input.value = ''; // Clear textboxes
                       }
                   });
               }
           });
       }

       function toggleInputs(section, disabled) {
           // Disable/enable all input fields within a section
           const inputs = section.querySelectorAll("input");
           inputs.forEach(input => {
               input.disabled = disabled;
           });
       }    
     
       function toggleSoftware(selectedSoftware) {
           const softwareIds = ['nwchem', 'qchem', 'psi4'];
       
           softwareIds.forEach(software => {
               const checkbox = document.getElementById(`software_${software}`);
               const optionsDiv = document.getElementById(`${software}-options`);
       
               if (software === selectedSoftware) {
                   // Toggle visibility of selected software options
                   if (checkbox.checked) {
                       optionsDiv.style.display = "block";
                   } else {
                       optionsDiv.style.display = "none";
                   }
               } else {
                   // Hide other software options and uncheck their boxes
                   const otherCheckbox = document.getElementById(`software_${software}`);
                   otherCheckbox.checked = false;
                   optionsDiv.style.display = "none";
               }
           });
       }
       function toggleOptionsEnabled(optionsDiv, enabled) {
           const inputs = optionsDiv.querySelectorAll("input, select, textarea");
           inputs.forEach(input => {
               input.disabled = !enabled;
           });
       }
       




       function generateFile(event) {
           // Gather form inputs
           event.preventDefault();

           const filename = document.getElementById("filename").value;
           const pdb_file = document.getElementById("pdb_file").value;
           const template_path = document.getElementById("template_path").value;
           const cutoff_radius = document.getElementById("cutoff_radius").value;
           const seed_ligand = document.getElementById("seed_ligand").checked;
           const seed_id = document.getElementById("seed_id").value;
           const seed_file = document.getElementById("seed_file").value;
           const charge_scheme = document.getElementById("charge_scheme").value;
           const software = document.getElementById("software_psi4").checked
               ? "psi4"
               : document.getElementById("software_nwchem").checked
               ? "nwchem"
               : document.getElementById("software_qchem").checked
               ? "q-chem"
               : null;
           const ligand_charge = document.getElementById("ligand_charge").value;
           const method = document.getElementById("method").value;
           const basis_set = document.getElementById("basis_set").value;
           const waterModel = document.getElementById("water_model").value;
           const oCharge = document.getElementById("o_charge").value;
           const hCharge = document.getElementById("h_charge").value;
           const epCharge = document.getElementById("ep_charge").value;
           const forcefield = document.getElementById("forcefield_charmm").checked
               ? "CHARMM"
               : document.getElementById("forcefield_amber").checked
               ? "Amber"
               : null;
           let capped = false;
           if (forcefield === "Amber") {
                capped = document.getElementById("precapbox").checked ? "false" : "true";
           }
           const charmm_rtf = document.getElementById("charmm_rtf").value;
           const charmm_prm = document.getElementById("charmm_prm").value;
           const cpChecked = document.getElementById("cp_correction").checked ? "true" : "false";
           const amber_ff = document.getElementById("amber_ff").value;
           const other_amber_ffs = document.getElementById("other_amber_ffs").value;
           const nwchem_scratch = document.getElementById("nwchem_scratch").value;
           const nwchem_perm = document.getElementById("nwchem_perm").value;
           const nwchem_scf = document.getElementById("nwchem_scf").value;
           const nwchem_dft = document.getElementById("nwchem_dft").value;
           const nwchem_mem = document.getElementById("mem").value;
           const qcsapt = document.getElementById("qchem_options").value;
           const qcopt = document.getElementById("qchem_sapt").value;
           const fisapt_partition = document.getElementById("fisapt_partition").checked ? "true" : "false";
           const do_fsapt = document.getElementById("do_fsapt").checked ? "true" : "false";
           const psi4options = document.getElementById("psi4options").value;
           const num_threads = document.getElementById("num_threads").value;
           const memory = document.getElementById("memory").value;
           // Validation logic
           if ((oCharge && !hCharge) || (!oCharge && hCharge)) {
               alert("If you enter either an Hydrogen or an Oxygen Charge, you must enter both.");
               return;
           }

           if (epCharge && (!oCharge || !hCharge)) {
               alert("If you enter an Extra Point Charge, you must enter both Oxygen Charge and Hydrogen Charge.");
               return;
           }


           const charmmRtf = document.getElementById("charmm_rtf").value;
           const charmmPrm = document.getElementById("charmm_prm").value;
           if (document.getElementById("forcefield_charmm").checked) {
               if (!charmmRtf || !charmmPrm) {
                   alert("Please fill in both CHARMM RTF and CHARMM PRM when CHARMM is selected.");
                   return;
               }
           }

           // Ensure Amber FF is filled if Amber is checked
           const amberFF = document.getElementById("amber_ff").value;
           if (document.getElementById("forcefield_amber").checked && !amberFF) {
               alert("Please fill in Amber FF when Amber is selected.");
               return;
           }


           // Validate Seed ID and Seed File when Ligand is unchecked
           if (!seed_ligand && (!seed_id || !seed_file)) {
               alert("Please fill in both seed_id and seed_file or check the Ligand box.");
               return;
           }

           // Create the content for the file
           let seed_content = seed_ligand
               ? "seed: ligand"
               : `seed_id: ${seed_id}\nseed_file: ${seed_file}`;

           let templateContent = template_path ? `\ntemplate_path: ${template_path}` : '';
           let hContent = hCharge ? `\nh_charge: ${hCharge}` : '';
           let oContent = oCharge ? `\no_charge: ${oCharge}` : '';
           let epContent = epCharge ? `\nep_charge: ${epCharge}` : '';
           let amber1 = amber_ff ? `\namber_ff: ${amber_ff}` : '';
           let amber2 = other_amber_ffs ? `\nother_amber_ffs: ${other_amber_ffs}` : '';
           let charmmp = charmm_prm ? `\ncharmm_prm: ${charmm_prm}` : '';
           let charmmr = charmm_rtf ? `\ncharmm_rtf: ${charmm_rtf}` : '';
           let nwcscratch = (software === 'nwchem' && nwchem_scratch) ? `\nnwchem_scratch: ${nwchem_scratch}` : '';
           let nwcperm = (software === 'nwchem' && nwchem_perm) ? `\nnwchem_perm: ${nwchem_perm}` : '';
           let nwcscf = (software === 'nwchem' && nwchem_scf) ? `\nnwchem_scf: ${nwchem_scf}` : '';
           let nwcdft = (software === 'nwchem' && nwchem_dft) ? `\nnwchem_dft: ${nwchem_dft}` : '';
           let nwcmem = (software === 'nwchem' && nwchem_mem) ? `\nnwchem_mem: ${nwchem_mem}` : '';
           let qopt = (software === 'q-chem' && qcopt) ? `\nqchem_options: ${qcopt}` : '';
           let qsapt = (software === 'q-chem' && qcsapt) ? `\nqchem_sapt: ${qcsapt}` : '';
           let fsapt_part = (software === 'psi4' && fisapt_partition && method === 'fisapt0') ? `\nfisapt_partition: ${fisapt_partition}` : '';
           let dofsapt = (software === 'psi4' && do_fsapt && method === 'fisapt0') ? `\ndo_fsapt: ${do_fsapt}` : '';
           let psi4opt = (software === 'psi4' && psi4options) ? `\npsi4_options: ${psi4options}` : '';
           let cp = (!method.includes('sapt') && cpChecked) ? `\ncp: ${cpChecked}` : ''; 
           let cap = capped ? `\npre-capped: ${capped}` : ''; 
           let nthreads = num_threads ? `\nnum_threads: ${num_threads}` : '';
           const content = `pdb_file: ${pdb_file}${templateContent}
   cutoff_radius: ${cutoff_radius}${amber1}${amber2}${cap}${charmmp}${charmmr}
   ${seed_content}
   charge_scheme: ${charge_scheme}
   water_model: ${waterModel}${hContent}${oContent}${epContent}
   software: ${software}${nwcscratch}${nwcperm}${nwcscf}${qopt}${qsapt}${fsapt_part}${dofsapt}${psi4opt}${nthreads}${cp}
   ligand_charge: ${ligand_charge}
   method: ${method}
   basis_set: ${basis_set}`;

           const amberChecked = document.getElementById('forcefield_amber').checked;
           const charmmChecked = document.getElementById('forcefield_charmm').checked;
       
           if (!amberChecked && !charmmChecked) {
               alert('Please select either Amber or CHARMM model type.');
               return; 
           }
       
           let message = '';
       
           if (charmmChecked) {
               message = 'Ensure the PDB is protein+environment, that the ligand is in the working directory as ligand.pdb, and there is a psf file of the protein in the same directory as the PDB';
               alert(message);
           }
           const blob = new Blob([content], { type: "text/plain" });
           const url = URL.createObjectURL(blob);

           const a = document.createElement("a");
           a.href = url;
           a.download = filename ? filename : "output.txt";
           document.body.appendChild(a);
           a.click();
           document.body.removeChild(a);

           URL.revokeObjectURL(url);
       }

       document.getElementById("seed_ligand").addEventListener("change", function() {
           const seedIdField = document.getElementById("seed_id");
           const seedFileField = document.getElementById("seed_file");

           if (this.checked) {
               // Disable the fields when Ligand checkbox is checked
               seedIdField.disabled = true;
               seedFileField.disabled = true;
               seedIdField.value = ''; // Clear seed_id field if ligand is checked
               seedFileField.value = ''; // Clear seed_file field if ligand is checked

               // Apply a darker background to show it's disabled
               seedIdField.style.backgroundColor = "#f0f0f0";
               seedFileField.style.backgroundColor = "#f0f0f0";
           } else {
               // Enable the fields when Ligand checkbox is unchecked
               seedIdField.disabled = false;
               seedFileField.disabled = false;

               // Reset background color to default
               seedIdField.style.backgroundColor = "";
               seedFileField.style.backgroundColor = "";
           }
       });
       document.getElementById("charmm-options").style.display = "none";
       document.getElementById("amber-options").style.display = "none";
   </script>

