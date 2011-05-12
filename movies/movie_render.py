from pymol import stored, cmd, selector, util
session = str(os.environ["SESSION"])
frames = int(os.environ["FRAMES"])
renderers = int(os.environ["RENDERERS"])
renderer_num = int(os.environ["PBS_ARRAYID"])
renderer_num = renderer_num - 1
print "STEP 1 open session."
cmd.load(session)
cmd.set("orthoscopic", 1)
cmd.set("ray_shadows", 1)
cmd.set("depth_cue", 1)
cmd.set("ray_trace_fog", 1)
cmd.set("antialias", 1.0)
cmd.set("cartoon_ring_mode", 3)
print "STEP 2 init movie"
cmd.mset("1 x%i" %(frames))
print "STEP 3 roll the movie."
#movie.roll(1, frames, 1, 'y')
util.mroll(1, frames, 1, 'y')
frames_per_renderer = frames / renderers
start_frame = 1 + (renderer_num * frames_per_renderer)
end_frame = start_frame + frames_per_renderer
print "TESTME FRAMES PER RENDERER: " + str(frames_per_renderer)
print "TESTME START FRAME THIS JOB: " + str(start_frame)
print "TESTME END FRAME THIS JOB: " + str(end_frame)
cmd.set("cache_frames", 0)
cmd.set("ray_trace_fog", 1)
cmd.set("ray_trace_frames", 1)
cmd.mpng(session, start_frame, end_frame)



