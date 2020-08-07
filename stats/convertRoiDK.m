function [out] = convertRoiDK(roi)

    out = '';
    switch (roi)
        case 'unknown'
            out = 'Unknown';
        case 'bankssts'
            out = 'Bank of Superior Temporal Sulcus';
        case 'caudalanteriorcingulate'
            out = 'Caudal Anterior Cingulate';
        case 'caudalmiddlefrontal'
            out = 'Caudal Middle Frontal';
        case 'corpuscallosum'
            out = 'Corpus Callosum';
        case 'cuneus'
            out = 'Cuneus';
        case 'entorhinal'
            out = 'Entorhinal';
        case 'fusiform'
            out = 'Fusiform';
        case 'inferiorparietal'
            out = 'Inferior Parietal';
        case 'inferiortemporal'
            out = 'Inferior Temporal';
        case 'isthmuscingulate'
            out = 'Isthmus Cingulate';
        case 'lateraloccipital'
            out = 'Lateral Occipital';
        case 'lateralorbitofrontal'
            out = 'Lateral Orbitofrontal';
        case 'lingual'
            out = 'Lingual';
        case 'medialorbitofrontal'
            out = 'Medial Orbitofrontal';
        case 'middletemporal'
            out = 'Middle Temporal';
        case 'parahippocampal'
            out = 'Parahippocampal';
        case 'paracentral'
            out = 'Paracentral';
        case 'parsopercularis'
            out = 'Pars Opercularis';
        case 'parsorbitalis'
            out = 'Pars Orbitalis';
        case 'parstriangularis'
            out = 'Pars Triangularis';
        case 'pericalcarine'
            out = 'Pericalcarine';
        case 'postcentral'
            out = 'Postcentral';
        case 'posteriorcingulate'
            out = 'Posterior Cingulate';
        case 'precentral'
            out = 'Precentral';
        case 'precuneus'
            out = 'Precuneus';
        case 'rostralanteriorcingulate'
            out = 'Rostral Anterior Cingulate';
        case 'rostralmiddlefrontal'
            out = 'Rostral Middle Frontal';
        case 'superiorfrontal'
            out = 'Superior Frontal';
        case 'superiorparietal'
            out = 'Superior Parietal';
        case 'superiortemporal'
            out = 'Superior Temporal';
        case 'supramarginal'
            out = 'Supramarginal';
        case 'frontalpole'
            out = 'Frontal Pole';
        case 'temporalpole'
            out = 'Temporal Pole';
        case 'transversetemporal'
            out = 'Transverse Temporal';
        case 'insula'
            out = 'Insula';
    end

end

