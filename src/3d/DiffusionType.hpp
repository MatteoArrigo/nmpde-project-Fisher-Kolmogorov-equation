#ifndef DIFFUSION_TYPE_H
#define DIFFUSION_TYPE_H


/**
 *
 * The DiffusionType enum is used to distinguish between different anisotropic diffusion behaviors:
 * - ISOTROPIC: Uniform diffusion in all directions.
 * - CIRCUMFERENTIAL: Diffusion primarily in the circumferential direction.
 * - RADIAL: Diffusion primarily in the radial direction.
 * 
 * Note: We did not implement the axon-based anisotropic diffusion model, due to lack of real-world data.
 */


enum DiffusionType {
    ISOTROPIC       = 0,
    CIRCUMFERENTIAL = 1,
    RADIAL          = 2
};

// Function to get the name of the diffusion type
inline std::string DiffusionName(DiffusionType type) {
    switch (type) {
        case ISOTROPIC:       return "ISOTROPIC";
        case CIRCUMFERENTIAL: return "CIRCUMFERENTIAL";
        case RADIAL:          return "RADIAL";
        default:              return "UNKNOWN";
    }
}

#endif // DIFFUSION_TYPE_H
