from paraview.simple import *
import glob
import os
import csv
import sys


# ============================================================
# FUNCTION
# ============================================================

def calculate_contour(folder, dim, data):

    # ========================================================
    # CHECK INPUTS
    # ========================================================

    if dim not in [2, 3]:
        raise ValueError("ERROR: dim must be either 2 or 3. Code stopped.")

    if not os.path.isdir(folder):
        raise ValueError(f"ERROR: Folder does not exist: {folder}")

    if not data:
        raise ValueError("ERROR: Data folder name cannot be empty.")

    
    # ========================================================
    # OUTPUT FOLDER
    # ========================================================

    output_folder = os.path.join(folder, data)

    os.makedirs(output_folder, exist_ok=True)

    print("Input folder :", folder)
    print("Output folder:", output_folder)
    print("Dimension    :", dim)


    # ========================================================
    # FIND FILES
    # ========================================================

    pvtu_pattern = os.path.join(folder, "*.pvtu")
    vtu_pattern = os.path.join(folder, "*.vtu")

    pvtu_files = sorted(glob.glob(pvtu_pattern))
    vtu_files = sorted(glob.glob(vtu_pattern))

    # ========================================================
    # SELECT FILE TYPE AND READER
    # ========================================================

    if len(pvtu_files) > 0:

        files = pvtu_files

        print("Detected partitioned VTU data.")
        print("Using .pvtu files.")
        print("Found", len(files), "files")

        solution = XMLPartitionedUnstructuredGridReader(
            registrationName="Solution", FileName=files
        )

    elif len(vtu_files) > 0:

        files = vtu_files

        print("Detected standard VTU data.")
        print("Using .vtu files.")
        print("Found", len(files), "files")

        solution = XMLUnstructuredGridReader(
            registrationName="Solution", FileName=files
        )

    else:

        raise RuntimeError(
            "ERROR: No .pvtu or .vtu files found!"
        )


    # ========================================================
    # OUTPUT FILE
    # ========================================================

    OUTPUT_FILE = os.path.join(
        output_folder,
        f"Contour{dim}D.csv"
    )


    # ========================================================
    # READ FILES
    # ========================================================

    solution.UpdatePipeline()

    # ========================================================
    # TIME STEPS
    # ========================================================

    animationScene = GetAnimationScene()
    animationScene.UpdateAnimationUsingDataTimeSteps()

    timesteps = animationScene.TimeKeeper.TimestepValues

    print("Number of timesteps:", len(timesteps))


    # ========================================================
    # CREATE CONTOUR
    # ========================================================

    contour1 = Contour(Input=solution)

    contour1.ContourBy = ["POINTS", "c"]

    contour1.Isosurfaces = [0.5]

    contour1.UpdatePipeline()


    # ========================================================
    # OPEN CSV
    # ========================================================

    if os.path.exists(OUTPUT_FILE):
        print("WARNING: Output file already exists and will be overwritten:")
        print(OUTPUT_FILE)

    with open(OUTPUT_FILE, "w", newline="") as csvfile:

        writer = csv.writer(csvfile)

        if dim == 2:
            writer.writerow([
                "TimeStep", "Time", "X", "Y", "Z",
                "c", "p", "u", "v",
                "velocity_x", "velocity_y", "velocity_z"
            ])

        elif dim == 3:
            writer.writerow([
                "TimeStep", "Time", "X", "Y", "Z",
                "c", "p", "u", "v", "w",
                "velocity_x", "velocity_y", "velocity_z"
            ])


        # ====================================================
        # LOOP OVER TIME
        # ====================================================

		data_found = False
		
        for i, time in enumerate(timesteps):

            #print("Processing timestep", i, "/", len(timesteps) - 1)

            solution.UpdatePipeline(time)
            contour1.UpdatePipeline(time)

            data_contour = servermanager.Fetch(contour1)

            if data_contour is None:
                print("WARNING: No data")
                continue

            number_of_points = data_contour.GetNumberOfPoints()

            #print("Contour points:", number_of_points)

            if number_of_points == 0:
                print("WARNING: Empty contour at timestep", i)
                continue
                
            data_found = True

            # =================================================
            # DATA ARRAYS
            # =================================================

            point_data = data_contour.GetPointData()

            c_array = point_data.GetArray("c")
            p_array = point_data.GetArray("p")
            u_array = point_data.GetArray("u")
            v_array = point_data.GetArray("v")

            if dim == 3:
                w_array = point_data.GetArray("w")

            velocity_array = point_data.GetArray("velocity")


            # =================================================
            # WRITE POINTS
            # =================================================

            for point_id in range(number_of_points):

                x, y, z = data_contour.GetPoint(point_id)


                # ---------------------------------------------
                # Scalars
                # ---------------------------------------------

                c = c_array.GetTuple1(point_id) if c_array is not None else ""

                p = p_array.GetTuple1(point_id) if p_array is not None else ""

                u = u_array.GetTuple1(point_id) if u_array is not None else ""

                v = v_array.GetTuple1(point_id) if v_array is not None else ""

                if dim == 3:
                    w = w_array.GetTuple1(point_id) if w_array is not None else ""


                # ---------------------------------------------
                # Velocity vector
                # ---------------------------------------------

                if velocity_array is not None:

                    velocity = velocity_array.GetTuple(point_id)

                    velocity_x = velocity[0]
                    velocity_y = velocity[1]

                    if len(velocity) > 2:
                        velocity_z = velocity[2]
                    else:
                        velocity_z = 0.0

                else:

                    velocity_x = ""
                    velocity_y = ""
                    velocity_z = ""


                # ---------------------------------------------
                # Write row
                # ---------------------------------------------

                if dim == 2:

                    writer.writerow([
                        i, time, x, y, z, c, p, u, v,
                        velocity_x, velocity_y, velocity_z
                    ])

                elif dim == 3:

                    writer.writerow([
                        i, time, x, y, z, c, p, u, v, w,
                        velocity_x, velocity_y, velocity_z
                    ])

	# ========================================================
    # CHECK DATA
    # ========================================================

    if not data_found:
        print("")
        print("WARNING: No contour data was found.")
        
    # ========================================================
    # CLEANUP
    # ========================================================

    Delete(contour1)
    Delete(solution)


    # ========================================================
    # FINISHED
    # ========================================================

    print("")
    print("============================================")
    if data_found:
        print("DONE")
    else:
        print("DONE - NO CONTOUR DATA")
    print("============================================")
    print("Files found :", len(files))
    print("Output file :", OUTPUT_FILE)
    print("============================================")

	if data_found:
        return OUTPUT_FILE
    else:
        return None


# ============================================================
# SCALABILITY
# ============================================================

def StrongScalability(folder, step):

    # ========================================================
    # CHECK INPUTS
    # ========================================================

    if not os.path.isdir(folder):
        raise ValueError(f"ERROR: Folder does not exist: {folder}")

    # ========================================================
    # FIND PARALLEL CASES
    # ========================================================

    parallel_pattern = os.path.join(
        folder,
        "Solution-quad-lev4_Parallel_*"
    )

    parallel_folders = sorted(glob.glob(parallel_pattern))

    if len(parallel_folders) == 0:
        raise RuntimeError(
            "ERROR: No parallel simulation folders were found!"
        )

    # ========================================================
    # CHECK FILES
    # ========================================================

    cases = []

    for parallel_folder in parallel_folders:

        if not os.path.isdir(parallel_folder):
            continue

        folder_name = os.path.basename(parallel_folder)

        ranks = int(folder_name.split("_")[-1])

        integral_file = os.path.join(
            parallel_folder,
            "Integral.txt"
        )

        if not os.path.isfile(integral_file):
            raise RuntimeError(
                f"ERROR: Integral.txt does not exist: {integral_file}"
            )

        # ====================================================
        # READ INTEGRAL
        # ====================================================

        numstep = None
        total_work_time = None

        with open(integral_file, "r") as file:

            for line in file:

                if line.startswith("Step") or line.startswith("-"):
                    continue

                columns = line.split()

                if len(columns) < 6:
                    continue

                current_step = int(columns[0])

                if current_step == step:

                    numstep = current_step
                    total_work_time = float(columns[5])

        if total_work_time is None:
            raise RuntimeError(
                f"ERROR: Step {step} was not found in: {integral_file}"
            )

        print(
            "Ranks:", ranks,
            " Step:", numstep,
            " Total work time:", total_work_time
        )

        cases.append([
            ranks,
            numstep,
            total_work_time,
            parallel_folder,
            integral_file
        ])

    # ========================================================
    # FINISHED
    # ========================================================

    print("")

    print("============================================")
    print("SCALABILITY FILE CHECK")
    print("============================================")
    print("Cases found:", len(cases))
    print("Step:", step)
    print("All Integral.txt files exist.")
    print("============================================")

    return cases
# ============================================================
# WEAK SCALABILITY
# ============================================================

def WeakScalability(folder, step, factor):

    # ========================================================
    # CHECK INPUTS
    # ========================================================

    if not os.path.isdir(folder):
        raise ValueError(f"ERROR: Folder does not exist: {folder}")

    if factor <= 0:
        raise ValueError("ERROR: factor must be greater than zero.")

    # ========================================================
    # FIND PARALLEL CASES
    # ========================================================

    parallel_pattern = os.path.join(
        folder,
        "Solution-quad-lev*_Parallel_*"
    )

    parallel_folders = sorted(glob.glob(parallel_pattern))

    if len(parallel_folders) == 0:
        raise RuntimeError(
            "ERROR: No parallel simulation folders were found!"
        )

    # ========================================================
    # LEVELS
    # ========================================================

    levels = [1, 2, 3, 5]

    # ========================================================
    # FIND INITIAL LEVEL AND n0
    # ========================================================

    lev0 = None
    n0 = None

    for parallel_folder in parallel_folders:

        if not os.path.isdir(parallel_folder):
            continue

        folder_name = os.path.basename(parallel_folder)

        lev = int(
            folder_name.split("_")[0].split("lev")[-1]
        )

        parallel_rank = int(
            folder_name.split("_")[-1]
        )

        if lev not in levels:
            continue

        if lev0 is None or lev < lev0:

            lev0 = lev
            n0 = parallel_rank

    if lev0 is None:
        raise RuntimeError(
            "ERROR: No valid weak scalability levels were found!"
        )

    # ========================================================
    # CHECK FILES
    # ========================================================

    cases = []

    for parallel_folder in parallel_folders:

        if not os.path.isdir(parallel_folder):
            continue

        folder_name = os.path.basename(parallel_folder)

        lev = int(
            folder_name.split("_")[0].split("lev")[-1]
        )

        if lev not in levels:
            continue

        parallel_rank = int(
            folder_name.split("_")[-1]
        )

        expected_rank = n0 * factor**(lev - lev0)

        if parallel_rank != expected_rank:
            continue

        integral_file = os.path.join(
            parallel_folder,
            "Integral.txt"
        )

        if not os.path.isfile(integral_file):
            raise RuntimeError(
                f"ERROR: Integral.txt does not exist: {integral_file}"
            )

        # ====================================================
        # READ INTEGRAL
        # ====================================================

        numstep = None
        total_work_time = None

        with open(integral_file, "r") as file:

            for line in file:

                if line.startswith("Step") or line.startswith("-"):
                    continue

                columns = line.split()

                if len(columns) < 6:
                    continue

                current_step = int(columns[0])

                if current_step == step:

                    numstep = current_step
                    total_work_time = float(columns[5])

        if total_work_time is None:
            raise RuntimeError(
                f"ERROR: Step {step} was not found in: {integral_file}"
            )

        print(
            "Lev:", lev,
            " Parallel ranks:", parallel_rank,
            " Step:", numstep,
            " Total work time:", total_work_time
        )

        cases.append([
            lev,
            parallel_rank,
            numstep,
            total_work_time,
            parallel_folder,
            integral_file
        ])

    # ========================================================
    # FINISHED
    # ========================================================

    print("")

    print("============================================")
    print("WEAK SCALABILITY FILE CHECK")
    print("============================================")
    print("Cases found:", len(cases))
    print("Step:", step)
    print("Factor:", factor)
    print("Initial level:", lev0)
    print("Initial ranks:", n0)
    print("============================================")

    return cases


# ============================================================
# COMMAND-LINE INTERFACE
# ============================================================

if __name__ == "__main__":

    if len(sys.argv) < 2:
        raise RuntimeError(
            "Usage: pvpython FunctionTools.py <function> <arguments>"
        )

    function = sys.argv[1]

    if function == "contour":

        if len(sys.argv) != 5:
            raise RuntimeError(
                "Usage: pvpython FunctionTools.py contour <folder> <dim> <data>"
            )

        folder = sys.argv[2]
        dim = int(sys.argv[3])
        data = sys.argv[4]

        calculate_contour(folder, dim, data)

    elif function == "strong_scalability":

        if len(sys.argv) != 4:
            raise RuntimeError(
                "Usage: pvpython FunctionTools.py strong scalability <folder> <step>"
            )

        folder = sys.argv[2]
        step = int(sys.argv[3])

        StrongScalability(folder, step)
    elif function == "weak_scalability":
        if len(sys.argv) != 5:
            raise RuntimeError("Usage: pvpython FunctionTools.py weak_scalability <folder> <step> <factor>")

        folder = sys.argv[2]
        step = int(sys.argv[3])
        factor = int(sys.argv[4])
        
        WeakScalability(folder, step, factor)

    else:

        raise RuntimeError(
            f"ERROR: Unknown function: {function}"
        )
