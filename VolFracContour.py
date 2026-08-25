from paraview.simple import *
import glob
import os
import csv

# ============================================================
# PATH
# ============================================================

folder = (
    "/Users/gonzald/Desktop/ECCOMAS/"
    "SolutionAvanche1Parallel2D_limex_L2H1L2/"
    "Solution_0000-quad-lev4-Rusanov-RelVel-NoForce-AveDiff-SlipDiff"
)

pattern = os.path.join(folder, "*.pvtu")

files = sorted(glob.glob(pattern))

OUTPUT_FILE = os.path.join(
    folder,
    "Contour2D.csv"
)


# ============================================================
# FIND FILES
# ============================================================

files = sorted(glob.glob(pattern))

print("Found", len(files), "files")

if len(files) == 0:
    raise RuntimeError("No .pvtu files found!")


# ============================================================
# READ FILES
# ============================================================

solution = XMLPartitionedUnstructuredGridReader(
    registrationName="Solution",
    FileName=files
)

solution.UpdatePipeline()


# ============================================================
# TIME STEPS
# ============================================================

animationScene = GetAnimationScene()
animationScene.UpdateAnimationUsingDataTimeSteps()

timesteps = animationScene.TimeKeeper.TimestepValues

print("Number of timesteps:", len(timesteps))


# ============================================================
# CREATE CONTOUR
# ============================================================

contour1 = Contour(Input=solution)

contour1.ContourBy = [
    "POINTS",
    "c"
]

contour1.Isosurfaces = [0.5]

contour1.UpdatePipeline()


# ============================================================
# OPEN CSV
# ============================================================

with open(
    OUTPUT_FILE,
    "w",
    newline=""
) as csvfile:

    writer = csv.writer(csvfile)

    writer.writerow([
        "TimeStep",
        "Time",
        "X",
        "Y",
        "Z",
        "c",
        "p",
        "u",
        "v",
        "velocity_x",
        "velocity_y",
        "velocity_z"
    ])


    # ========================================================
    # LOOP OVER TIME
    # ========================================================

    for i, time in enumerate(timesteps):

        print(
            "Processing timestep",
            i,
            "/",
            len(timesteps) - 1
        )

        solution.UpdatePipeline(time)
        contour1.UpdatePipeline(time)

        data = servermanager.Fetch(contour1)

        if data is None:
            print("WARNING: No data")
            continue

        number_of_points = data.GetNumberOfPoints()

        print(
            "Contour points:",
            number_of_points
        )

        if number_of_points == 0:
            print(
                "WARNING: Empty contour at timestep",
                i
            )
            continue


        # ====================================================
        # DATA ARRAYS
        # ====================================================

        point_data = data.GetPointData()

        c_array = point_data.GetArray("c")
        p_array = point_data.GetArray("p")
        u_array = point_data.GetArray("u")
        v_array = point_data.GetArray("v")
        velocity_array = point_data.GetArray("velocity")


        # ====================================================
        # WRITE POINTS
        # ====================================================

        for point_id in range(number_of_points):

            x, y, z = data.GetPoint(point_id)


            # -----------------------------------------------
            # Scalars
            # -----------------------------------------------

            c = (
                c_array.GetTuple1(point_id)
                if c_array is not None
                else ""
            )

            p = (
                p_array.GetTuple1(point_id)
                if p_array is not None
                else ""
            )

            u = (
                u_array.GetTuple1(point_id)
                if u_array is not None
                else ""
            )

            v = (
                v_array.GetTuple1(point_id)
                if v_array is not None
                else ""
            )


            # -----------------------------------------------
            # Velocity vector
            # -----------------------------------------------

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


            # -----------------------------------------------
            # Write row
            # -----------------------------------------------

            writer.writerow([
                i,
                time,
                x,
                y,
                z,
                c,
                p,
                u,
                v,
                velocity_x,
                velocity_y,
                velocity_z
            ])


# ============================================================
# FINISHED
# ============================================================

print("")
print("============================================")
print("DONE")
print("============================================")
print("Files found :", len(files))
print("Output file :", OUTPUT_FILE)
print("============================================")