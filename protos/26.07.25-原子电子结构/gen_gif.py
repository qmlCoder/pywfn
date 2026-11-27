from PIL import Image

paths=[f'./imgs/elems_{i+1:0>3}.png' for i in range(118)]

images=[Image.open(path) for path in paths]

images[0].save('elems.gif', save_all=True, append_images=images[1:], duration=100, loop=0)
