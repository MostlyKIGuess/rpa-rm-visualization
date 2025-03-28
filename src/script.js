const m = 3;
const r = 2;
const n = 1 << m;
let Nmax = 1;
let initialY = [0, 0, 1, 0, 0, 0, 0, 0];

const rm1_m2_codewords = [
  [0, 0, 0, 0],
  [1, 1, 1, 1],
  [0, 1, 0, 1],
  [1, 0, 1, 0],
  [0, 0, 1, 1],
  [1, 1, 0, 0],
  [0, 1, 1, 0],
  [1, 0, 0, 1],
].map((arr) => new Uint8Array(arr));

const ziIndices = Array.from({ length: n - 1 }, (_, i) => i + 1);

const nmaxInput = document.getElementById("nmaxInput");
const startButton = document.getElementById("startButton");
const stepsContainer = document.getElementById("steps");
const vectorInputs = document.querySelectorAll(".bit-input");
const currentVectorDisplay = document.getElementById("currentVectorDisplay");

function updateButtonText() {
  startButton.textContent = `Run Visualization (Nₘₐₓ=${Nmax})`;
}

// Update Nmax when input changes
nmaxInput.addEventListener("input", function () {
  Nmax = parseInt(this.value) || 1;
  if (Nmax < 1) Nmax = 1;
  if (Nmax > 10) Nmax = 10;
  this.value = Nmax;
  updateButtonText();
});

function updateVectorDisplay() {
  const vector = Array.from(vectorInputs).map(
    (input) => parseInt(input.value) || 0
  );
  initialY = vector;
  currentVectorDisplay.textContent = vector.join("");
}

function updateVectorDisplay() {
  const vector = Array.from(vectorInputs).map(
    (input) => parseInt(input.value) || 0
  );
  initialY = vector;
  currentVectorDisplay.textContent = vector.join("");
}

vectorInputs.forEach((input) => {
  input.addEventListener("input", function () {
    if (this.value !== "0" && this.value !== "1") {
      this.value = "0";
    }
    updateVectorDisplay();
  });
});

function intToBinVec(val, bits) {
  const binStr = val.toString(2).padStart(bits, "0");
  return new Uint8Array(binStr.split("").map(Number));
}

function binVecToInt(vec) {
  return parseInt(vec.join(""), 2);
}

function vecToStr(vec) {
  return vec.join("");
}

function hammingDistance(v1, v2) {
  let dist = 0;
  for (let i = 0; i < v1.length; i++) {
    if (v1[i] !== v2[i]) {
      dist++;
    }
  }
  return dist;
}

function getCosetMapping(ziVec, m_local) {
  const n_local = 1 << m_local;
  const n_proj = 1 << (m_local - 1);
  const zi_idx = binVecToInt(ziVec);

  let processed = new Uint8Array(n_local);
  let coset_map = {};
  let coset_pairs = {};

  let complement_basis_indices = [];
  if (m_local !== 3) {
    let basis_count = 0;
    for (let i = 0; i < m_local && basis_count < m_local - 1; ++i) {
      const basis_idx = 1 << i;
      if (basis_idx !== zi_idx) {
        complement_basis_indices.push(basis_idx);
        basis_count++;
      }
    }
    if (basis_count !== m_local - 1) return null;
  } else {
    if (zi_idx === 1) complement_basis_indices = [4, 2];
    else if (zi_idx === 2) complement_basis_indices = [4, 1];
    else if (zi_idx === 4) complement_basis_indices = [2, 1];
    else if (zi_idx === 3) complement_basis_indices = [4, 1];
    else if (zi_idx === 5) complement_basis_indices = [2, 1];
    else if (zi_idx === 6) complement_basis_indices = [2, 1];
    else if (zi_idx === 7) complement_basis_indices = [4, 2];
    else return null;
  }

  let proj_idx_k = 0;
  for (let ab = 0; ab < 1 << (m_local - 1); ++ab) {
    const ab_vec = intToBinVec(ab, m_local - 1);
    let rep_idx = 0;
    for (let bit_idx = 0; bit_idx < m_local - 1; ++bit_idx) {
      if (ab_vec[bit_idx] === 1) {
        rep_idx ^= complement_basis_indices[bit_idx];
      }
    }

    if (!processed[rep_idx]) {
      const idx1 = rep_idx;
      const idx2 = rep_idx ^ zi_idx;

      coset_map[idx1] = ab;
      coset_pairs[ab] = [idx1, idx2];

      processed[idx1] = 1;
      processed[idx2] = 1;
    }
  }

  if (processed.reduce((a, b) => a + b, 0) !== n_local) {
    let unprocessed_indices = [];
    for (let i = 0; i < n_local; ++i) {
      if (!processed[i]) unprocessed_indices.push(i);
    }
    return null;
  }

  return { coset_map, coset_pairs };
}

function project(y, zi_idx, m_local, addLogEntry) {
  const n_local = 1 << m_local;
  const n_proj = 1 << (m_local - 1);
  const zi_vec = intToBinVec(zi_idx, m_local);
  let proj_vector = new Uint8Array(n_proj);

  addLogEntry(
    `<h4>Projecting onto \\(B_{${zi_idx}} = \\{0, ${vecToStr(zi_vec)}\\}\\)</h4>`,
    false
  );
  const mapping = getCosetMapping(zi_vec, m_local);
  if (!mapping) {
    addLogEntry(
      `<p class="explanation">Error generating coset mapping for zi=${vecToStr(
        zi_vec
      )}. Skipping projection.</p>`,
      true
    );
    return null;
  }
  const { coset_pairs } = mapping;

  let logDetails =
    '<p class="explanation">Calculating XOR for each coset index k (0 to 3):</p><ul>';
  for (let k = 0; k < n_proj; k++) {
    if (!coset_pairs[k]) {
      logDetails += `<li>Error: No pair found for k=${k}</li>`;
      continue;
    }
    const [idx1, idx2] = coset_pairs[k];
    const val1 = y[idx1];
    const val2 = y[idx2];
    const xor_sum = val1 ^ val2;
    proj_vector[k] = xor_sum;
    logDetails += `<li>k=${k} (orig idx ${idx1}='${vecToStr(
      intToBinVec(idx1, m_local)
    )}', ${idx2}='${vecToStr(
      intToBinVec(idx2, m_local)
    )}'): y[${idx1}] \\(\\oplus\\) y[${idx2}] = ${val1} \\(\\oplus\\) ${val2} = ${xor_sum}</li>`;
  }
  logDetails += "</ul>";
  addLogEntry(logDetails, false);

  addLogEntry(
    `<p>Result <span class="vector">y/B<sub>${zi_idx}</sub></span> = <span class="vector">${vecToStr(
      proj_vector
    )}</span></p>`,
    true
  );
  return proj_vector;
}

function fhtDecodeRm1(y_proj, m_proj, rm1_codewords_local, addLogEntry) {
  const n_proj_local = 1 << m_proj;
  let min_dist = n_proj_local + 1;
  let decoded_codeword = null;
  let ties = [];

  addLogEntry(
    `<h4>FHT Decoding RM(${m_proj},1) for <span class="vector">${vecToStr(
      y_proj
    )}</span></h4>`,
    false
  );

  let distanceTable =
    "<div class='table-container'><p>Hamming distances to all codewords:</p><table>";
  distanceTable += "<tr><th>Codeword</th><th>Hamming Distance</th></tr>";

  for (const codeword of rm1_codewords_local) {
    const dist = hammingDistance(y_proj, codeword);

    let rowClass = "";
    if (dist < min_dist) {
      min_dist = dist;
      decoded_codeword = codeword;
      ties = [codeword];
      rowClass = "highlight-row";
    } else if (dist === min_dist) {
      ties.push(codeword);
      rowClass = "highlight-row";
    }

    distanceTable += `<tr class="${rowClass}"><td>${vecToStr(
      codeword
    )}</td><td>${dist}</td></tr>`;
  }
  distanceTable += "</table></div>";
  addLogEntry(distanceTable, false);

  if (ties.length > 1) {
    addLogEntry(
      `<p class="explanation">Note: Found ${
        ties.length
      } closest codewords (dist=${min_dist}). Using first found: ${vecToStr(
        decoded_codeword
      )}.</p>`,
      false
    );
  } else {
    addLogEntry(
      `<p class="explanation">Unique closest codeword found (dist=${min_dist}).</p>`,
      false
    );
  }

  addLogEntry(
    `<p>Decoded <span class="vector">ŷ/B</span> = <span class="vector">${vecToStr(
      decoded_codeword
    )}</span></p>`,
    true
  );
  return new Uint8Array(decoded_codeword);
}

function aggregate(y, y_hat_proj_map, m_local, addLogEntry) {
  const n_local = 1 << m_local;
  const n_proj = 1 << (m_local - 1);
  let changevote = new Int32Array(n_local);

  addLogEntry("<h3>Step 3: Aggregation (Algorithm 2)</h3>", false);
  addLogEntry(
    `<p class="explanation">Combine results using majority vote. Input y = <span class="vector">${vecToStr(
      y
    )}</span>.</p>`,
    false
  );

  for (let z_idx = 0; z_idx < n_local; z_idx++) {
    for (let zi_idx = 1; zi_idx < n_local; zi_idx++) {
      const zi_vec = intToBinVec(zi_idx, m_local);
      const mapping = getCosetMapping(zi_vec, m_local);
      if (!mapping) continue;
      const { coset_pairs } = mapping;

      let proj_idx_k = -1;
      for (let k_test = 0; k_test < n_proj; k_test++) {
        if (!coset_pairs[k_test]) continue;
        if (
          z_idx === coset_pairs[k_test][0] ||
          z_idx === coset_pairs[k_test][1]
        ) {
          proj_idx_k = k_test;
          break;
        }
      }
      if (proj_idx_k === -1) continue;

      const [idx1, idx2] = coset_pairs[proj_idx_k];
      const y_proj_val = y[idx1] ^ y[idx2];
      const decoded_proj = y_hat_proj_map[zi_idx];

      if (!decoded_proj) continue;
      const y_hat_proj_val = decoded_proj[proj_idx_k];

      if (y_proj_val !== y_hat_proj_val) {
        changevote[z_idx]++;
      }
    }
  }

  let tableHTML =
    "<h4>Calculated Change Votes:</h4><table><tr><th>z (Index)</th><th>z (Binary)</th><th>changevote(z)</th></tr>";
  for (let z_idx = 0; z_idx < n_local; z_idx++) {
    const isOverThreshold = changevote[z_idx] > (n_local - 1) / 2.0;
    tableHTML += `<tr ${
      isOverThreshold ? 'class="highlight-row"' : ""
    }><td>${z_idx}</td><td>${vecToStr(intToBinVec(z_idx, m_local))}</td><td>${
      changevote[z_idx]
    }</td></tr>`;
  }
  tableHTML += "</table>";
  addLogEntry(tableHTML, false);

  let y_new = new Uint8Array(y);
  const threshold = (n_local - 1) / 2.0;
  addLogEntry(
    `<p class="explanation">Majority Threshold = (n-1)/2 = ${
      n_local - 1
    }/2 = ${threshold.toFixed(
      1
    )}. Flip bit z if changevote(z) > ${threshold.toFixed(
      1
    )} (i.e., \(\geq 4\)).</p>`,
    false
  );

  let flipped_indices = [];
  for (let z_idx = 0; z_idx < n_local; z_idx++) {
    if (changevote[z_idx] > threshold) {
      y_new[z_idx] = y[z_idx] ^ 1;
      flipped_indices.push(z_idx);
    }
  }

  if (flipped_indices.length > 0) {
    addLogEntry(
      `<p class="explanation">Flipping bits at original indices: ${flipped_indices.join(
        ", "
      )}.</p>`,
      false
    );
  } else {
    addLogEntry(
      `<p class="explanation">No changevote count exceeded the threshold. No bits flipped.</p>`,
      false
    );
  }

  addLogEntry(
    `<p>New estimate <span class="vector">ŷ</span> = <span class="vector">${vecToStr(
      y_new
    )}</span></p>`,
    true
  );
  return { y_new, changevote };
}

function runVisualization() {
  stepsContainer.innerHTML = "";
  MathJax.typeset();

  let currentY = new Uint8Array(initialY);
  let visualizationLog = [];
  let iterationCount = 0;
  let converged = false;

  const addLogEntry = (html, isMajorStep) => {
    if (isMajorStep) {
      const stepDiv = document.createElement("div");
      stepDiv.className = "step";
      stepDiv.innerHTML = visualizationLog.join("");
      stepsContainer.appendChild(stepDiv);
      visualizationLog = [html];
      MathJax.typeset([stepDiv]);
    } else {
      visualizationLog.push(html);
    }
  };

  while (iterationCount < Nmax && !converged) {
    iterationCount++;

    addLogEntry(`<h2>Iteration 1 / ${Nmax}</h2>`, true);
    addLogEntry(
      `<p>Current vector estimate <span class="vector">y</span> = <span class="vector">${vecToStr(
        currentY
      )}</span></p>`,
      false
    );

    let y_proj_map = {};
    let y_hat_proj_map = {};

    addLogEntry(
      "<h3>Step 1 & 2: Projection & Recursive Decode (Base Case FHT)</h3>",
      true
    );
    addLogEntry(
      `<p class="explanation">For each of the \\(n-1 = 7\\) non-zero vectors \\(z_i\\), project \\(y\\) onto the cosets of \\(B_i=\\{0, z_i\\}\\) to get \\(y/B_i\\). Then, decode \\(y/B_i\\) using FHT for \\(RM(m-1=2, r-1=1)\\).</p>`,
      false
    );

    for (const zi_idx of ziIndices) {
      const proj_vector = project(currentY, zi_idx, m, addLogEntry);
      if (proj_vector === null) {
        addLogEntry(
          `<p class="explanation">Skipping decoding for \(B_{${zi_idx}}\) due to projection error.</p>`,
          true
        );
        continue;
      }
      y_proj_map[zi_idx] = proj_vector;

      const decoded_proj = fhtDecodeRm1(
        proj_vector,
        m - 1,
        rm1_m2_codewords,
        addLogEntry
      );
      y_hat_proj_map[zi_idx] = decoded_proj;
    }

    const { y_new, changevote } = aggregate(
      currentY,
      y_hat_proj_map,
      m,
      addLogEntry
    );

    addLogEntry("<h3>Step 4: Convergence Check</h3>", true);
    converged = y_new.every((val, index) => val === currentY[index]);

    addLogEntry(
      `<p class="explanation">Compare y before aggregation (<span class="vector">${vecToStr(
        currentY
      )}</span>) with ŷ after aggregation (<span class="vector">${vecToStr(
        y_new
      )}</span>).</p>`,
      false
    );

    if (converged) {
      addLogEntry(
        '<p class="explanation highlight">Result: \\(y = \\hat{y}\\). Algorithm Converged!</p>',
        false
      );
    } else {
      addLogEntry(
        '<p class="explanation highlight">Result: \\(y \\neq \\hat{y}\\). Algorithm Did Not Converge in this iteration.</p>',
        false
      );
    }

    const finalStepDiv = document.createElement("div");
    finalStepDiv.className = "step";
    finalStepDiv.innerHTML = visualizationLog.join("");
    stepsContainer.appendChild(finalStepDiv);
    visualizationLog = [];

    MathJax.typeset([finalStepDiv]);

    addLogEntry("<h2>Final Result</h2>", true);
    addLogEntry(
      `<p class="explanation">After ${Nmax} iteration(s), the algorithm ${
        converged ? "converged" : "stopped"
      }.</p>`,
      false
    );
    addLogEntry(
      `<p>Final decoded codeword estimate <span class="vector">ĉ</span> = <span class="vector">${vecToStr(
        y_new
      )}</span></p>`,
      false
    );

    // Compare with actual transmitted codeword
    const transmittedCodeword = new Uint8Array(n).fill(0); // All zeros
    const success = y_new.every((val, idx) => val === transmittedCodeword[idx]);

    if (success) {
      addLogEntry(
        `<p class="explanation highlight">Success! The decoded codeword matches the transmitted codeword \\(c=(00000000)\\).</p>`,
        false
      );
    } else {
      addLogEntry(
        `<p class="explanation highlight">The decoding failed. The decoded codeword \\(\\hat{c}=(${vecToStr(
          y_new
        )})\\) does not match the transmitted codeword \\(c=(00000000)\\).</p>`,
        false
      );
    }

    const finalResultDiv = document.createElement("div");
    finalResultDiv.className = "step";
    finalResultDiv.innerHTML = visualizationLog.join("");
    stepsContainer.appendChild(finalResultDiv);

    MathJax.typeset([finalResultDiv]);
    if (converged) break;
  }
  // Scroll to results
  stepsContainer.scrollIntoView({ behavior: "smooth", block: "start" });
}

startButton.addEventListener("click", runVisualization);

// Initialize the vector display
updateVectorDisplay();
updateButtonText();

// Make sure MathJax typesets the initial content
document.addEventListener("DOMContentLoaded", function () {
  setTimeout(() => {
    MathJax.typeset();
  }, 1000);
});
