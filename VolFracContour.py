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

    pattern = os.path.join(folder, "*.pvtu")

    files = sorted(glob.glob(pattern))

    print("Found", len(files), "files")

    if len(files) == 0:
        raise RuntimeError("No .pvtu files found!")


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

    solution = XMLPartitionedUnstructuredGridReader(
        registrationName="Solution", FileName=files
    )

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

        for i, time in enumerate(timesteps):

            print("Processing timestep", i, "/", len(timesteps) - 1)

            solution.UpdatePipeline(time)
            contour1.UpdatePipeline(time)

            data_contour = servermanager.Fetch(contour1)

            if data_contour is None:
                print("WARNING: No data")
                continue

            number_of_points = data_contour.GetNumberOfPoints()

            print("Contour points:", number_of_points)

            if number_of_points == 0:
                print("WARNING: Empty contour at timestep", i)
                continue


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
    # CLEANUP
    # ========================================================

    Delete(contour1)
    Delete(solution)


    # ========================================================
    # FINISHED
    # ========================================================

    print("")
    print("============================================")
    print("DONE")
    print("============================================")
    print("Files found :", len(files))
    print("Output file :", OUTPUT_FILE)
    print("============================================")

    return OUTPUT_FILE


# ============================================================
# COMMAND-LINE INTERFACE
# ============================================================

if __name__ == "__main__":

    if len(sys.argv) != 4:
        raise RuntimeError(
            "Usage: pvpython VolFracContour.py <folder> <dim> <data>"
        )

    folder = sys.argv[1]
    dim = int(sys.argv[2])
    data = sys.argv[3]

    calculate_contour(folder, dim, data)
