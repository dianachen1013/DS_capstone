% Define the path to the Schaefer 2018 parcellation NIfTI file
schaefer_nifti_path = 'D:/ibp/schaefer/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii';

% Load left and right hemispheres
[surf_lh, surf_rh] = load_conte69();

% Load the Schaefer 2018 atlas (read the NIfTI file)
schaefer_nifti = load_nii(schaefer_nifti_path);
schaefer_atlas = Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii;

% Verify the atlas dimensions
disp(size(schaefer_atlas)); % Should match your brain volume dimensions

% Map the atlas to the cortical surface (parcellate it)
labeling = schaefer_atlas_to_surface(surf_lh, surf_rh, schaefer_atlas);

% Define clusters (replace these with your actual ROI indices)
cluster1 = [1 36 38 39 49 68 69 77 85 86 87 88 89 90 92 94 95 96 102 103 104 112 113 114 115 116 117 118 119 120 121 122 123 124 126 127 131 132 134 135 137 139 140 141 147 148 149 150 153 154 157 158 159 161 162 165 166 167 169 171 174 177 178 180 181 182 184 187 195 232 238 242 271 274 290 295 296 297 299 300 303 304 305 308 309 319 321 322 323 326 327 328 330 331 332 333 337 338 339 340 342 343 344 346 347 349 350 353 354 359 362 363 364 365 366 367 368 369 370 371 372 373 374 375 376 377 378 380 382 383 385 386 387 388 389 390 393 399];
cluster2 = [2 3 4 5 7 8 10 11 13 14 16 17 18 22 23 26 27 28 30 41 51 52 56 58 60 70 71 72 73 74 75 78 79 80 81 93 138 145 160 168 170 172 175 179 183 201 203 204 205 206 208 210 211 213 215 218 221 222 223 226 227 229 275 280 282 284];
cluster3 = [33 34 35 37 40 42 43 44 45 46 47 48 50 53 54 55 57 59 61 62 63 64 65 66 67 83 84 91 100 101 105 111 136 142 151 152 155 156 188 209 216 228 230 235 236 237 239 240 241 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 272 273 276 278 279 283 286 287 288 289 293 298 301 306 307 310 312 314 315 316 317 318 335 341 348 351 352];
cluster4 = [6 9 12 15 19 20 21 24 25 29 31 32 76 82 97 98 99 106 107 108 109 110 125 128 129 130 133 143 144 146 163 164 173 176 185 186 189 190 191 192 193 194 196 197 198 199 200 202 207 212 214 217 219 220 224 225 231 233 234 277 281 285 291 292 294 302 311 313 320 324 325 329 334 336 345 355 356 357 358 360 361 379 381 384 391 392 394 395 396 397 398 400];

% Initialize a label vector for the 400 ROIs
cluster_labels = zeros(400, 1); % 400 ROIs in Schaefer atlas

% Assign cluster indices to the label vector
cluster_labels(cluster1) = 1; % Label cluster1 as 1
cluster_labels(cluster2) = 2; % Label cluster2 as 2
cluster_labels(cluster3) = 3; % Label cluster3 as 3
cluster_labels(cluster4) = 4; % Label cluster4 as 4

% Plot your clusters on the cortical surface
plot_hemispheres(cluster_labels, {surf_lh, surf_rh}, ...
    'parcellation', labeling, ...
    'overlay_type', 'label', ...
    'colormap', [1 0 0; 0 1 0; 0 0 1; 1 1 0]);  % Example colormap for 4 clusters: Red, Green, Blue, Yellow

title('Brain Surface Visualization of 4 Clusters');
view([-90 0]);  % Adjust the view as needed
colorbar;
