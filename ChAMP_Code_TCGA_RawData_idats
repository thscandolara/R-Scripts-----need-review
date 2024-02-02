#código para aplicar o champ nos dados brutos (.idats) derivados de amostras de laringe do TCGA;

> myLoad <- champ.load(directory = getwd(),
+            method = "ChAMP",
+           arraytype = "450k",
+           methValue = "B",
+           autoimpute = F,
+           filterDetP = T,
+           ProbeCutoff = 0,
+           SampleCutoff = 0.1,
+           detPcut = 0.01,
+           filterBeads = T,
+           beadCutoff = 0.05,
+           filterNoCG = F,
+           filterSNPs = F,
+           population = NULL,
+           filterMultiHit = T,
+           filterXY = F,
+           force = F)
[===========================]
[<<<< ChAMP.LOAD START >>>>>]
-----------------------------

[ Loading Data with ChAMP Method ]
----------------------------------
Note that ChAMP method will NOT return rgSet or mset, they object defined by minfi. Which means, if you use ChAMP method to load data, you can not use SWAN or FunctionNormliazation method in champ.norm() (you can use BMIQ or PBC still). But All other function should not be influenced.

[===========================]
[<<<< ChAMP.IMPORT START >>>>>]
-----------------------------

[ Section 1: Read PD Files Start ]
  CSV Directory: /data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/sample_sheet.csv
  Find CSV Success
  Reading CSV File
  Replace Sentrix_Position into Array
  Replace Sentrix_ID into Slide
  There is NO Pool_ID in your pd file.
  There is NO Sample_Plate in your pd file.
  There is NO Sample_Well in your pd file.
[ Section 1: Read PD file Done ]


[ Section 2: Read IDAT files Start ]
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/507d6335-4dc5-4ff2-b75a-768b50fcbe64_noid_Grn.idat ---- (1/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/427bfdcc-0562-420c-82c9-fad413ec534f_noid_Grn.idat ---- (2/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/a675ad13-5314-407e-8869-5669351f4140_noid_Grn.idat ---- (3/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/d05b00bf-afbd-492c-a12c-decb30ef1c21_noid_Grn.idat ---- (4/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/d6bcb720-f239-4fda-ae38-82fd84e9fbac_noid_Grn.idat ---- (5/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/cd749253-c3c0-4b7f-942a-dedf37d69f39_noid_Grn.idat ---- (6/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/63530c8b-9719-4417-8619-84503e6929c3_noid_Grn.idat ---- (7/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/878330e3-1cfd-44f3-8110-4ac477bf96b7_noid_Grn.idat ---- (8/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/71543d6c-ddc1-42d5-8ef8-c1a92bc56c9e_noid_Grn.idat ---- (9/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/841411a2-c9f1-4259-999b-bb32aa0ea3e0_noid_Grn.idat ---- (10/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/9fccbe69-d3fe-4211-bee9-27f52eb8f00f_noid_Grn.idat ---- (11/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/0685483a-e5cc-42f0-a09f-b03cf0dc875a_noid_Grn.idat ---- (12/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/0cc5492e-a8c5-4571-bfe4-818d7617d649_noid_Grn.idat ---- (13/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/a283ec51-4868-49b3-b48f-9e717aac46dc_noid_Grn.idat ---- (14/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/5941e79f-101f-4bc2-99f4-7017c60a068f_noid_Grn.idat ---- (15/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/4435d721-ef8d-4224-8b06-2d0ef22e8530_noid_Grn.idat ---- (16/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/5bc63ec4-849b-4176-900a-5beb5e86a020_noid_Grn.idat ---- (17/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/dbdff281-0a2e-49bf-bc48-7b980d944e69_noid_Grn.idat ---- (18/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ef3655ee-909e-4b87-a195-eaf7f6e85be5_noid_Grn.idat ---- (19/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/43e1fe57-0614-4b31-bf20-0dcec95229b3_noid_Grn.idat ---- (20/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/f2ce72ee-e3b4-418f-93ee-62a8f719b74e_noid_Grn.idat ---- (21/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/0550903a-7e6a-4351-98a6-42cdabc292dc_noid_Grn.idat ---- (22/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/4b3445b1-0a5c-4844-b7b6-d7847877d548_noid_Grn.idat ---- (23/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/7ad846cd-e494-4d64-8fa5-a9af1418a6c4_noid_Grn.idat ---- (24/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/6ffe162a-def0-4f83-9691-f169aa0f081a_noid_Grn.idat ---- (25/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/116b08d3-b565-4b1f-aee3-643cfbdf4d3d_noid_Grn.idat ---- (26/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/7155880f-003a-46db-82ef-17f3025fd451_noid_Grn.idat ---- (27/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/fbd0a9cf-ac46-4d7f-9b31-eab9bb0304dc_noid_Grn.idat ---- (28/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/2d1646e5-f0c5-4ca0-a4ef-990113c53c1c_noid_Grn.idat ---- (29/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/67ed70fb-1cef-4e05-8b9e-d50344d0640e_noid_Grn.idat ---- (30/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ebbfbb58-4d91-4b1b-b0f1-ffa7183b24bb_noid_Grn.idat ---- (31/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/a4b331a9-ef16-4fde-a881-8ec06d69f48c_noid_Grn.idat ---- (32/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/a5251926-a56a-4cc7-9e39-1786bf182afa_noid_Grn.idat ---- (33/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/8d5b9023-2250-4df9-a1fd-f25831768544_noid_Grn.idat ---- (34/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ff93e589-ca51-43f3-b988-b2bfacf6b4b9_noid_Grn.idat ---- (35/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/d8cbd220-fea3-4e55-b799-3a082a5e0ae8_noid_Grn.idat ---- (36/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/64e6337a-b8bb-4ff2-9bf4-729378bf1dae_noid_Grn.idat ---- (37/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/d71a5700-1e65-4014-8903-157fc8312545_noid_Grn.idat ---- (38/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/f120c4ee-f8da-4b8a-8bf8-cb2d461494c1_noid_Grn.idat ---- (39/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/f74e9ae2-77ec-4304-8a98-b4dd24756100_noid_Grn.idat ---- (40/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/09c33348-5341-4482-b6fc-fe99fb1984f6_noid_Grn.idat ---- (41/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/2e3312a2-6f6b-47eb-92b1-39cf43ecad78_noid_Grn.idat ---- (42/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/94236e6e-5788-4435-a8b4-3cc1ed77f3db_noid_Grn.idat ---- (43/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ef678599-ff02-4c08-a5ff-1277eb93974d_noid_Grn.idat ---- (44/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/e1645f10-7fe5-4c67-83d0-d85930f1fcde_noid_Grn.idat ---- (45/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ab97ae87-ce0f-4a1c-a89f-f980b8b7f6e5_noid_Grn.idat ---- (46/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/1df48486-67b3-400f-a1be-1f64f7c9ead6_noid_Grn.idat ---- (47/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/583c30ac-f52c-4994-9c6b-9aa78fc3c500_noid_Grn.idat ---- (48/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/3c99dae5-c443-44f6-a606-a0e9f56e08f9_noid_Grn.idat ---- (49/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/80f2dfb6-d3ea-42be-8f02-4a0035652fc0_noid_Grn.idat ---- (50/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/b9d196c9-a302-4c26-ad72-69c929f6cc18_noid_Grn.idat ---- (51/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/a42b961a-7691-4dee-ab6a-894dc411ca41_noid_Grn.idat ---- (52/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/c9fef126-0e97-4301-ad90-e1fdedf34acb_noid_Grn.idat ---- (53/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/776119c3-4546-4f3e-91e8-eb18e3b27244_noid_Grn.idat ---- (54/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ff6a91f4-11ef-4d78-8412-143db9669e85_noid_Grn.idat ---- (55/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/3afeb448-ad90-4333-896c-9a022239c281_noid_Grn.idat ---- (56/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/9ad943a1-3de9-4a09-866c-4510b5aa0f90_noid_Grn.idat ---- (57/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/77a4e246-7566-4a0d-8172-70d1af4d7d01_noid_Grn.idat ---- (58/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/06496da4-2928-4714-a2d8-dfc5418d5d0a_noid_Grn.idat ---- (59/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/aad9d915-fe1b-4fb3-b31a-91c4dc2f77b8_noid_Grn.idat ---- (60/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ca6763a4-6293-4b8f-bf32-b7208fe3be5b_noid_Grn.idat ---- (61/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/40a24e46-b32f-41ce-8890-9690dff93f52_noid_Grn.idat ---- (62/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/733af45e-272c-4a7b-8440-18a5cc736dad_noid_Grn.idat ---- (63/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/0646ebf0-8afa-4f77-adb8-18ddf17fa367_noid_Grn.idat ---- (64/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/bd1d9066-c8c6-4a4f-8c6b-e0784995ab01_noid_Grn.idat ---- (65/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/546ab9c1-bd0d-4a7b-ab49-25645de867f9_noid_Grn.idat ---- (66/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/a91245a8-81c2-4a01-813f-30f7a3d64c31_noid_Grn.idat ---- (67/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/0bf5e63e-edcb-44d6-8aba-0ceec9b586ae_noid_Grn.idat ---- (68/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/7dec7522-3763-4eff-bc64-04bae89f82d9_noid_Grn.idat ---- (69/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/267c5e01-7fdb-48fc-912f-0216f16e194f_noid_Grn.idat ---- (70/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/0ec5e19e-86a3-47fb-8929-8d32a389a89a_noid_Grn.idat ---- (71/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/809b43ce-d4be-4df5-a0a0-a2227beed544_noid_Grn.idat ---- (72/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ed4f6ea3-4fa6-4325-b2ab-3a77fce60649_noid_Grn.idat ---- (73/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/20f2bbda-9a43-49aa-a5c4-8b4d7ecabdbe_noid_Grn.idat ---- (74/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/c0af9291-33f9-471b-b344-2a93f02654aa_noid_Grn.idat ---- (75/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/fa815826-5b90-4d72-ad9a-451245e6f79d_noid_Grn.idat ---- (76/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ec87298d-6ce6-41b1-a639-26d9d0cdca44_noid_Grn.idat ---- (77/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/758045b1-3dea-40c4-a442-7db227257a4d_noid_Grn.idat ---- (78/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/3e0e2943-e7ef-4305-9439-48c34bb5d07d_noid_Grn.idat ---- (79/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/2dffc286-5d6f-477f-86f3-48ff5ac51b1f_noid_Grn.idat ---- (80/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/490e45e5-5ec2-4caf-af20-425bceb172c1_noid_Grn.idat ---- (81/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/a7e4e9ed-7867-4a70-a0ca-d18eb294990c_noid_Grn.idat ---- (82/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/de8f47b7-5e76-463c-a252-b080e34ba872_noid_Grn.idat ---- (83/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/f002b61a-a5f7-43ac-b02c-eb0a629eba6e_noid_Grn.idat ---- (84/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/3322bd8f-4221-4178-b3ae-86d541af98bd_noid_Grn.idat ---- (85/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/fe9703e0-cefd-4c26-a750-90bebc6073dc_noid_Grn.idat ---- (86/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/921931bb-3386-40a5-808c-aa324d957022_noid_Grn.idat ---- (87/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/58f2cfa3-8dfd-4490-8872-ba7dcc6c8fd8_noid_Grn.idat ---- (88/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/c8e10e82-122d-4198-945c-5673a9959f2c_noid_Grn.idat ---- (89/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/b1a54bd1-d0ce-4f56-a638-830ec5e7a3f0_noid_Grn.idat ---- (90/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/eaa2be34-6db5-4f06-a60f-863e00b14db0_noid_Grn.idat ---- (91/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/93e4cc5f-b374-4c6b-ba2a-4522bbd84f2f_noid_Grn.idat ---- (92/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/e17ad668-733f-4e42-be6c-fc338bab31f8_noid_Grn.idat ---- (93/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/6e1e8e2c-e99e-4f1d-a97e-d1f034b05812_noid_Grn.idat ---- (94/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/0a0fe9ff-a0ad-4c1e-a3b8-ae9a0c3b3925_noid_Grn.idat ---- (95/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/38154d38-4c36-462b-ba4f-78fac5088a46_noid_Grn.idat ---- (96/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/97a3f792-c1aa-47be-8dac-98028030f00b_noid_Grn.idat ---- (97/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/9849e3b8-fefd-4aab-97c8-821bbce2b774_noid_Grn.idat ---- (98/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/55204b84-d175-4112-9eff-3e15916994c8_noid_Grn.idat ---- (99/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/b145357c-da52-4ab2-bf15-e58f1e5149df_noid_Grn.idat ---- (100/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/fb67f865-91bb-442a-aded-c7c6ffc629f3_noid_Grn.idat ---- (101/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/c3249ef1-d554-4fa9-8b7c-e7800ddc9b11_noid_Grn.idat ---- (102/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/11033594-1305-41c7-83ce-2cb702641a89_noid_Grn.idat ---- (103/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/5819df5d-0c3d-4821-b312-a7c3f3c42438_noid_Grn.idat ---- (104/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/0b459b25-3561-4cb5-a501-36b9ff8f729e_noid_Grn.idat ---- (105/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/4a825128-7489-45d8-9edd-cffc2b3733bb_noid_Grn.idat ---- (106/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/45ef47e2-7432-45be-9d0d-4df1f43d41f5_noid_Grn.idat ---- (107/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/aaf5e197-4360-46eb-9ab2-95b971783805_noid_Grn.idat ---- (108/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/471f0f73-1ad9-4f2b-8f62-dbcdd67d0ff1_noid_Grn.idat ---- (109/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/afd7e648-bc14-46c6-8b12-f25a234a8f4b_noid_Grn.idat ---- (110/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/5aafa689-51a4-458d-a642-3f1b3266741b_noid_Grn.idat ---- (111/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/f5318d95-34f1-4647-939a-3c9b47555c73_noid_Grn.idat ---- (112/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/af8ae291-7cb0-4675-a4c0-a4145d13f61b_noid_Grn.idat ---- (113/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/9bfc3487-c330-4c99-9cf6-b7bdfdd006ee_noid_Grn.idat ---- (114/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/d68501c5-973f-4d4a-baa4-9b347ee055e1_noid_Grn.idat ---- (115/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/4e0c878c-8d55-4f03-8251-ebab2a29c25e_noid_Grn.idat ---- (116/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/b4ca3291-69de-4494-a4f3-6b0c23bc7538_noid_Grn.idat ---- (117/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ace562fd-5d94-466b-ae3f-bc5e1e46ddd5_noid_Grn.idat ---- (118/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/c9bb9a1b-7499-4fda-8802-d6246d1b4a43_noid_Grn.idat ---- (119/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/02c81587-f4dc-4897-aa31-8ffdb9cfab2b_noid_Grn.idat ---- (120/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/8f4f1327-85c7-4a57-8e50-c836a0dc6a2b_noid_Grn.idat ---- (121/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/08d65ed4-749d-4ada-be8b-bde730f7b57f_noid_Grn.idat ---- (122/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/327643d8-231b-4b75-aa07-58ddd8fbbdfa_noid_Grn.idat ---- (123/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/3df99627-30c3-48f1-8148-2f9a28b27fab_noid_Grn.idat ---- (124/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ea533077-9d1d-480c-ac59-04021de8483a_noid_Grn.idat ---- (125/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/6010aeb4-61ce-4880-b710-459c6a5165a1_noid_Grn.idat ---- (126/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/78f8ceca-25bb-4549-b2c0-b0a65dda6b9c_noid_Grn.idat ---- (127/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/d7b80871-c027-4541-8276-94d116368ebf_noid_Grn.idat ---- (128/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/2297caa3-4bb0-4cba-8964-7bde7a65d155_noid_Grn.idat ---- (129/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/78ba68fd-204b-4aa0-94c5-246190c4dc2e_noid_Grn.idat ---- (130/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/507d6335-4dc5-4ff2-b75a-768b50fcbe64_noid_Red.idat ---- (1/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/427bfdcc-0562-420c-82c9-fad413ec534f_noid_Red.idat ---- (2/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/a675ad13-5314-407e-8869-5669351f4140_noid_Red.idat ---- (3/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/d05b00bf-afbd-492c-a12c-decb30ef1c21_noid_Red.idat ---- (4/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/d6bcb720-f239-4fda-ae38-82fd84e9fbac_noid_Red.idat ---- (5/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/cd749253-c3c0-4b7f-942a-dedf37d69f39_noid_Red.idat ---- (6/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/63530c8b-9719-4417-8619-84503e6929c3_noid_Red.idat ---- (7/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/878330e3-1cfd-44f3-8110-4ac477bf96b7_noid_Red.idat ---- (8/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/71543d6c-ddc1-42d5-8ef8-c1a92bc56c9e_noid_Red.idat ---- (9/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/841411a2-c9f1-4259-999b-bb32aa0ea3e0_noid_Red.idat ---- (10/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/9fccbe69-d3fe-4211-bee9-27f52eb8f00f_noid_Red.idat ---- (11/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/0685483a-e5cc-42f0-a09f-b03cf0dc875a_noid_Red.idat ---- (12/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/0cc5492e-a8c5-4571-bfe4-818d7617d649_noid_Red.idat ---- (13/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/a283ec51-4868-49b3-b48f-9e717aac46dc_noid_Red.idat ---- (14/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/5941e79f-101f-4bc2-99f4-7017c60a068f_noid_Red.idat ---- (15/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/4435d721-ef8d-4224-8b06-2d0ef22e8530_noid_Red.idat ---- (16/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/5bc63ec4-849b-4176-900a-5beb5e86a020_noid_Red.idat ---- (17/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/dbdff281-0a2e-49bf-bc48-7b980d944e69_noid_Red.idat ---- (18/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ef3655ee-909e-4b87-a195-eaf7f6e85be5_noid_Red.idat ---- (19/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/43e1fe57-0614-4b31-bf20-0dcec95229b3_noid_Red.idat ---- (20/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/f2ce72ee-e3b4-418f-93ee-62a8f719b74e_noid_Red.idat ---- (21/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/0550903a-7e6a-4351-98a6-42cdabc292dc_noid_Red.idat ---- (22/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/4b3445b1-0a5c-4844-b7b6-d7847877d548_noid_Red.idat ---- (23/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/7ad846cd-e494-4d64-8fa5-a9af1418a6c4_noid_Red.idat ---- (24/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/6ffe162a-def0-4f83-9691-f169aa0f081a_noid_Red.idat ---- (25/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/116b08d3-b565-4b1f-aee3-643cfbdf4d3d_noid_Red.idat ---- (26/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/7155880f-003a-46db-82ef-17f3025fd451_noid_Red.idat ---- (27/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/fbd0a9cf-ac46-4d7f-9b31-eab9bb0304dc_noid_Red.idat ---- (28/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/2d1646e5-f0c5-4ca0-a4ef-990113c53c1c_noid_Red.idat ---- (29/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/67ed70fb-1cef-4e05-8b9e-d50344d0640e_noid_Red.idat ---- (30/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ebbfbb58-4d91-4b1b-b0f1-ffa7183b24bb_noid_Red.idat ---- (31/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/a4b331a9-ef16-4fde-a881-8ec06d69f48c_noid_Red.idat ---- (32/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/a5251926-a56a-4cc7-9e39-1786bf182afa_noid_Red.idat ---- (33/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/8d5b9023-2250-4df9-a1fd-f25831768544_noid_Red.idat ---- (34/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ff93e589-ca51-43f3-b988-b2bfacf6b4b9_noid_Red.idat ---- (35/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/d8cbd220-fea3-4e55-b799-3a082a5e0ae8_noid_Red.idat ---- (36/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/64e6337a-b8bb-4ff2-9bf4-729378bf1dae_noid_Red.idat ---- (37/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/d71a5700-1e65-4014-8903-157fc8312545_noid_Red.idat ---- (38/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/f120c4ee-f8da-4b8a-8bf8-cb2d461494c1_noid_Red.idat ---- (39/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/f74e9ae2-77ec-4304-8a98-b4dd24756100_noid_Red.idat ---- (40/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/09c33348-5341-4482-b6fc-fe99fb1984f6_noid_Red.idat ---- (41/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/2e3312a2-6f6b-47eb-92b1-39cf43ecad78_noid_Red.idat ---- (42/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/94236e6e-5788-4435-a8b4-3cc1ed77f3db_noid_Red.idat ---- (43/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ef678599-ff02-4c08-a5ff-1277eb93974d_noid_Red.idat ---- (44/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/e1645f10-7fe5-4c67-83d0-d85930f1fcde_noid_Red.idat ---- (45/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ab97ae87-ce0f-4a1c-a89f-f980b8b7f6e5_noid_Red.idat ---- (46/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/1df48486-67b3-400f-a1be-1f64f7c9ead6_noid_Red.idat ---- (47/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/583c30ac-f52c-4994-9c6b-9aa78fc3c500_noid_Red.idat ---- (48/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/3c99dae5-c443-44f6-a606-a0e9f56e08f9_noid_Red.idat ---- (49/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/80f2dfb6-d3ea-42be-8f02-4a0035652fc0_noid_Red.idat ---- (50/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/b9d196c9-a302-4c26-ad72-69c929f6cc18_noid_Red.idat ---- (51/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/a42b961a-7691-4dee-ab6a-894dc411ca41_noid_Red.idat ---- (52/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/c9fef126-0e97-4301-ad90-e1fdedf34acb_noid_Red.idat ---- (53/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/776119c3-4546-4f3e-91e8-eb18e3b27244_noid_Red.idat ---- (54/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ff6a91f4-11ef-4d78-8412-143db9669e85_noid_Red.idat ---- (55/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/3afeb448-ad90-4333-896c-9a022239c281_noid_Red.idat ---- (56/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/9ad943a1-3de9-4a09-866c-4510b5aa0f90_noid_Red.idat ---- (57/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/77a4e246-7566-4a0d-8172-70d1af4d7d01_noid_Red.idat ---- (58/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/06496da4-2928-4714-a2d8-dfc5418d5d0a_noid_Red.idat ---- (59/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/aad9d915-fe1b-4fb3-b31a-91c4dc2f77b8_noid_Red.idat ---- (60/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ca6763a4-6293-4b8f-bf32-b7208fe3be5b_noid_Red.idat ---- (61/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/40a24e46-b32f-41ce-8890-9690dff93f52_noid_Red.idat ---- (62/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/733af45e-272c-4a7b-8440-18a5cc736dad_noid_Red.idat ---- (63/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/0646ebf0-8afa-4f77-adb8-18ddf17fa367_noid_Red.idat ---- (64/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/bd1d9066-c8c6-4a4f-8c6b-e0784995ab01_noid_Red.idat ---- (65/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/546ab9c1-bd0d-4a7b-ab49-25645de867f9_noid_Red.idat ---- (66/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/a91245a8-81c2-4a01-813f-30f7a3d64c31_noid_Red.idat ---- (67/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/0bf5e63e-edcb-44d6-8aba-0ceec9b586ae_noid_Red.idat ---- (68/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/7dec7522-3763-4eff-bc64-04bae89f82d9_noid_Red.idat ---- (69/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/267c5e01-7fdb-48fc-912f-0216f16e194f_noid_Red.idat ---- (70/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/0ec5e19e-86a3-47fb-8929-8d32a389a89a_noid_Red.idat ---- (71/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/809b43ce-d4be-4df5-a0a0-a2227beed544_noid_Red.idat ---- (72/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ed4f6ea3-4fa6-4325-b2ab-3a77fce60649_noid_Red.idat ---- (73/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/20f2bbda-9a43-49aa-a5c4-8b4d7ecabdbe_noid_Red.idat ---- (74/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/c0af9291-33f9-471b-b344-2a93f02654aa_noid_Red.idat ---- (75/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/fa815826-5b90-4d72-ad9a-451245e6f79d_noid_Red.idat ---- (76/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ec87298d-6ce6-41b1-a639-26d9d0cdca44_noid_Red.idat ---- (77/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/758045b1-3dea-40c4-a442-7db227257a4d_noid_Red.idat ---- (78/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/3e0e2943-e7ef-4305-9439-48c34bb5d07d_noid_Red.idat ---- (79/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/2dffc286-5d6f-477f-86f3-48ff5ac51b1f_noid_Red.idat ---- (80/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/490e45e5-5ec2-4caf-af20-425bceb172c1_noid_Red.idat ---- (81/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/a7e4e9ed-7867-4a70-a0ca-d18eb294990c_noid_Red.idat ---- (82/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/de8f47b7-5e76-463c-a252-b080e34ba872_noid_Red.idat ---- (83/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/f002b61a-a5f7-43ac-b02c-eb0a629eba6e_noid_Red.idat ---- (84/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/3322bd8f-4221-4178-b3ae-86d541af98bd_noid_Red.idat ---- (85/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/fe9703e0-cefd-4c26-a750-90bebc6073dc_noid_Red.idat ---- (86/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/921931bb-3386-40a5-808c-aa324d957022_noid_Red.idat ---- (87/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/58f2cfa3-8dfd-4490-8872-ba7dcc6c8fd8_noid_Red.idat ---- (88/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/c8e10e82-122d-4198-945c-5673a9959f2c_noid_Red.idat ---- (89/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/b1a54bd1-d0ce-4f56-a638-830ec5e7a3f0_noid_Red.idat ---- (90/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/eaa2be34-6db5-4f06-a60f-863e00b14db0_noid_Red.idat ---- (91/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/93e4cc5f-b374-4c6b-ba2a-4522bbd84f2f_noid_Red.idat ---- (92/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/e17ad668-733f-4e42-be6c-fc338bab31f8_noid_Red.idat ---- (93/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/6e1e8e2c-e99e-4f1d-a97e-d1f034b05812_noid_Red.idat ---- (94/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/0a0fe9ff-a0ad-4c1e-a3b8-ae9a0c3b3925_noid_Red.idat ---- (95/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/38154d38-4c36-462b-ba4f-78fac5088a46_noid_Red.idat ---- (96/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/97a3f792-c1aa-47be-8dac-98028030f00b_noid_Red.idat ---- (97/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/9849e3b8-fefd-4aab-97c8-821bbce2b774_noid_Red.idat ---- (98/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/55204b84-d175-4112-9eff-3e15916994c8_noid_Red.idat ---- (99/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/b145357c-da52-4ab2-bf15-e58f1e5149df_noid_Red.idat ---- (100/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/fb67f865-91bb-442a-aded-c7c6ffc629f3_noid_Red.idat ---- (101/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/c3249ef1-d554-4fa9-8b7c-e7800ddc9b11_noid_Red.idat ---- (102/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/11033594-1305-41c7-83ce-2cb702641a89_noid_Red.idat ---- (103/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/5819df5d-0c3d-4821-b312-a7c3f3c42438_noid_Red.idat ---- (104/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/0b459b25-3561-4cb5-a501-36b9ff8f729e_noid_Red.idat ---- (105/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/4a825128-7489-45d8-9edd-cffc2b3733bb_noid_Red.idat ---- (106/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/45ef47e2-7432-45be-9d0d-4df1f43d41f5_noid_Red.idat ---- (107/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/aaf5e197-4360-46eb-9ab2-95b971783805_noid_Red.idat ---- (108/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/471f0f73-1ad9-4f2b-8f62-dbcdd67d0ff1_noid_Red.idat ---- (109/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/afd7e648-bc14-46c6-8b12-f25a234a8f4b_noid_Red.idat ---- (110/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/5aafa689-51a4-458d-a642-3f1b3266741b_noid_Red.idat ---- (111/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/f5318d95-34f1-4647-939a-3c9b47555c73_noid_Red.idat ---- (112/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/af8ae291-7cb0-4675-a4c0-a4145d13f61b_noid_Red.idat ---- (113/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/9bfc3487-c330-4c99-9cf6-b7bdfdd006ee_noid_Red.idat ---- (114/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/d68501c5-973f-4d4a-baa4-9b347ee055e1_noid_Red.idat ---- (115/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/4e0c878c-8d55-4f03-8251-ebab2a29c25e_noid_Red.idat ---- (116/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/b4ca3291-69de-4494-a4f3-6b0c23bc7538_noid_Red.idat ---- (117/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ace562fd-5d94-466b-ae3f-bc5e1e46ddd5_noid_Red.idat ---- (118/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/c9bb9a1b-7499-4fda-8802-d6246d1b4a43_noid_Red.idat ---- (119/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/02c81587-f4dc-4897-aa31-8ffdb9cfab2b_noid_Red.idat ---- (120/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/8f4f1327-85c7-4a57-8e50-c836a0dc6a2b_noid_Red.idat ---- (121/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/08d65ed4-749d-4ada-be8b-bde730f7b57f_noid_Red.idat ---- (122/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/327643d8-231b-4b75-aa07-58ddd8fbbdfa_noid_Red.idat ---- (123/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/3df99627-30c3-48f1-8148-2f9a28b27fab_noid_Red.idat ---- (124/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/ea533077-9d1d-480c-ac59-04021de8483a_noid_Red.idat ---- (125/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/6010aeb4-61ce-4880-b710-459c6a5165a1_noid_Red.idat ---- (126/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/78f8ceca-25bb-4549-b2c0-b0a65dda6b9c_noid_Red.idat ---- (127/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/d7b80871-c027-4541-8276-94d116368ebf_noid_Red.idat ---- (128/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/2297caa3-4bb0-4cba-8964-7bde7a65d155_noid_Red.idat ---- (129/130)
  Loading:/data04/projects04/SheilaCoelho/headspace_methyl/data/idats_tcga_larynx_tumors/78ba68fd-204b-4aa0-94c5-246190c4dc2e_noid_Red.idat ---- (130/130)

  Extract Mean value for Green and Red Channel Success
    Your Red Green Channel contains 622399 probes.
[ Section 2: Read IDAT Files Done ]


[ Section 3: Use Annotation Start ]

  Reading 450k Annotation >>

  Fetching NEGATIVE ControlProbe.
    Totally, there are 613 control probes in Annotation.
    Your data set contains 613 control probes.

  Generating Meth and UnMeth Matrix
    Extracting Meth Matrix...
      Totally there are 485512 Meth probes in 450k Annotation.
      Your data set contains 485512 Meth probes.
    Extracting UnMeth Matrix...
      Totally there are 485512 UnMeth probes in 450k Annotation.
      Your data set contains 485512 UnMeth probes.

  Generating beta Matrix
  Generating M Matrix
  Generating intensity Matrix
  Calculating Detect P value
  Counting Beads
[ Section 3: Use Annotation Done ]

[<<<<< ChAMP.IMPORT END >>>>>>]
[===========================]
[You may want to process champ.filter() next.]

[===========================]
[<<<< ChAMP.FILTER START >>>>>]
-----------------------------

In New version ChAMP, champ.filter() function has been set to do filtering on the result of champ.import(). You can use champ.import() + champ.filter() to do Data Loading, or set "method" parameter in champ.load() as "ChAMP" to get the same effect.

This function is provided for user need to do filtering on some beta (or M) matrix, which contained most filtering system in champ.load except beadcount. User need to input beta matrix, pd file themselves. If you want to do filterintg on detP matrix and Bead Count, you also need to input a detected P matrix and Bead Count information.

Note that if you want to filter more data matrix, say beta, M, intensity... please make sure they have exactly the same rownames and colnames.


[ Section 1:  Check Input Start ]
  You have inputed beta,intensity for Analysis.

  pd file provided, checking if it's in accord with Data Matrix...
    pd file check success.

  Parameter filterDetP is TRUE, checking if detP in accord with Data Matrix...
    detP check success.

  Parameter filterBeads is TRUE, checking if beadcount in accord with Data Matrix...
    beadcount check success.

  Checking Finished :filterDetP,filterBeads,filterMultiHit would be done on beta,intensity.
  You also provided :detP,beadcount .
[ Section 1: Check Input Done ]


[ Section 2: Filtering Start >>

  Filtering Detect P value Start
    The fraction of failed positions per sample
    You may need to delete samples with high proportion of failed probes:

                 Failed CpG Fraction.
TCGA-F7-7848-01A         0.0007353062
TCGA-CN-4738-01A         0.0009659905
TCGA-CV-6962-01A         0.0187410404
TCGA-CV-7177-01A         0.0011266457
TCGA-CV-A460-01A         0.0046672379
TCGA-D6-6517-01A         0.0004345928
TCGA-CN-4735-01A         0.0012090329
TCGA-BB-4217-01A         0.0005210994
TCGA-D6-8568-01A         0.0003769217
TCGA-CN-A641-01A         0.0011843168
TCGA-CN-A49B-01A         0.0003007135
TCGA-F7-A50I-01A         0.0008053354
TCGA-CN-6012-01A         0.0027620327
TCGA-CN-4723-01A         0.0006570383
TCGA-CN-A63T-01A         0.0016209692
TCGA-CN-6010-01A         0.0010936908
TCGA-CV-5435-01A         0.0013758671
TCGA-CV-5441-01A         0.0002698183
TCGA-CV-5432-11B         0.0010421988
TCGA-T3-A92M-01A         0.0008918420
TCGA-CR-7399-01A         0.0014397172
TCGA-D6-A74Q-01A         0.0014088220
TCGA-CR-7374-01A         0.0012213910
TCGA-DQ-7595-01A         0.0003398474
TCGA-CV-5443-01A         0.0014232398
TCGA-CV-A45W-01A         0.0095486826
TCGA-CN-6988-01A         0.0087639440
TCGA-CV-5444-01A         0.0004963832
TCGA-CN-5360-01A         0.0021503073
TCGA-CV-7247-01A         0.0005499349
TCGA-UF-A7JJ-01A         0.0010648552
TCGA-CN-A497-01A         0.0002409827
TCGA-UF-A7J9-01A         0.0004222347
TCGA-BA-4076-01A         0.0020246667
TCGA-CR-7402-01A         0.0004757864
TCGA-CV-6962-11A         0.0193898400
TCGA-CR-6474-01A         0.0004181153
TCGA-CV-5430-01A         0.0024366030
TCGA-CV-5430-11B         0.0017115952
TCGA-QK-A8Z8-01A         0.0015571191
TCGA-CR-7370-01A         0.0011245860
TCGA-CN-6997-01A         0.0008794839
TCGA-CV-7433-01A         0.0004222347
TCGA-BA-A6DA-01A         0.0008032757
TCGA-CR-7388-01A         0.0008794839
TCGA-D6-6824-01A         0.0005561140
TCGA-CN-5361-01A         0.0009392147
TCGA-UF-A718-01A         0.0003542652
TCGA-UF-A7JK-01A         0.0004593089
TCGA-BA-5555-01A         0.0003007135
TCGA-CR-7364-01A         0.0004593089
TCGA-CV-7101-01A         0.0007867983
TCGA-CV-7101-11A         0.0005581736
TCGA-BB-7864-01A         0.0001915504
TCGA-DQ-5629-01A         0.0002203859
TCGA-CN-6989-01A         0.0142550545
TCGA-CV-7250-01A         0.0006549787
TCGA-CV-7250-11A         0.0003192506
TCGA-CV-6935-01A         0.0067022030
TCGA-CV-6935-11A         0.0104796586
TCGA-H7-A6C5-01A         0.0003913395
TCGA-CN-6023-01A         0.0027991069
TCGA-CV-7245-11A         0.0006487996
TCGA-BA-4078-01A         0.0028959119
TCGA-D6-A6EQ-01A         0.0003048328
TCGA-CN-4739-01A         0.0010319003
TCGA-CV-7261-01A         0.0003233700
TCGA-TN-A7HJ-01A         0.0011307650
TCGA-CV-7437-01A         0.0005231591
TCGA-BB-7862-01A         0.0003151312
TCGA-CN-A63W-01A         0.0021297105
TCGA-CN-6021-01A         0.0069102308
TCGA-CV-7410-01A         0.0003892798
TCGA-CR-7389-01A         0.0012131523
TCGA-CV-7089-01A         0.0007229481
TCGA-F7-A622-01A         0.0017939824
TCGA-D6-A6EK-01A         0.0002162665
TCGA-F7-8298-01A         0.0003377877
TCGA-CN-5355-01A         0.0019072649
TCGA-CV-7245-01A         0.0005684720
TCGA-BA-6868-01B         0.0081728155
TCGA-CN-A6V3-01A         0.0012996589
TCGA-CV-7440-01A         0.0005293381
TCGA-BA-6870-01A         0.0010133632
TCGA-HD-7229-01A         0.0004160556
TCGA-CN-4727-01A         0.0030194928
TCGA-DQ-7589-01A         0.0002904151
TCGA-CN-6992-01A         0.0057712271
TCGA-BA-6869-01A         0.0002306843
TCGA-UF-A71D-01A         0.0027888085
TCGA-CV-5443-11A         0.0027373165
TCGA-CV-7418-01A         0.0006076060
TCGA-BB-7870-01A         0.0004448912
TCGA-CV-5440-11A         0.0003089522
TCGA-CV-5440-01A         0.0004242944
TCGA-CN-A63U-01A         0.0003439668
TCGA-CN-5363-01A         0.0023665738
TCGA-D6-6826-01A         0.0074931207
TCGA-CV-7422-01A         0.0025643033
TCGA-CV-7248-01A         0.0012234507
TCGA-QK-A8ZB-01A         0.0005561140
TCGA-UF-A7JF-01A         0.0003274893
TCGA-CV-A6K1-01A         0.0002348037
TCGA-CV-5444-11A         0.0004654880
TCGA-CV-7242-01A         0.0016333273
TCGA-D6-A6ES-01A         0.0003213103
TCGA-CN-4722-01A         0.0022388736
TCGA-BA-A6DI-01A         0.0005437559
TCGA-CV-5435-11B         0.0006611577
TCGA-CV-5441-11A         0.0011884361
TCGA-CV-5432-01A         0.0015241642
TCGA-CV-5434-11B         0.0012090329
TCGA-CV-5434-01A         0.0012502266
TCGA-CV-7421-01A         0.0004860848
TCGA-CV-7430-01A         0.0006858739
TCGA-CV-5978-01A         0.0036003230
TCGA-CV-5978-11A         0.0036065020
TCGA-CR-7398-01A         0.0013676284
TCGA-CV-A45Z-01A         0.0047640429
TCGA-CV-7415-01A         0.0013243751
TCGA-CR-7371-01A         0.0010730940
TCGA-UF-A7JH-01A         0.0003357281
TCGA-CV-7089-11A         0.0009042001
TCGA-CV-A461-01A         0.0121377021
TCGA-F7-A623-01A         0.0016106708
TCGA-CN-6022-01A         0.0006055463
TCGA-CV-5431-11A         0.0006508593
TCGA-CV-5431-01A         0.0005911285
TCGA-CV-7424-01A         0.0009371550
TCGA-CN-5356-01A         0.0012646443

    Filtering probes with a detection p-value above 0.01.
    Removing 24361 probes.
    If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples

  Filtering BeadCount Start
    Filtering probes with a beadcount <3 in at least 5% of samples.
    Removing 592 probes

  Filtering MultiHit Start
    Filtering probes that align to multiple locations as identified in Nordlund et al
    Removing 8786 probes from the analysis.

  Updating PD file

  Fixing Outliers Start
    Replacing all value smaller/equal to 0 with smallest positive value.
    Replacing all value greater/equal to 1 with largest value below 1..
[ Section 2: Filtering Done ]

 All filterings are Done, now you have 451773 probes and 130 samples.

[<<<<< ChAMP.FILTER END >>>>>>]
[===========================]
[You may want to process champ.QC() next.]

[<<<<< ChAMP.LOAD END >>>>>>]
[===========================]
[You may want to process champ.QC() next.]
