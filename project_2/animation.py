import pygame
import sys

# Initialize Pygame
pygame.init()

# Screen dimensions
width, height = 1200, 800
screen = pygame.display.set_mode((width, height))

# Load background image
background_image = pygame.image.load(r"C:\Users\nullcline\school\spring_2024\bmeg422\bmeg422-project1\project_2\4.png").convert()
background_image = pygame.transform.scale(background_image, (width, height))

# Colors
colors = [  # Four different colors for the sweeps
    (255, 0, 0, 255),   # Red
    (0, 255, 0, 255),   # Green
    (0, 0, 255, 255),   # Blue
    (255, 255, 0, 255),  # Yellow
    (255,255, 255,)
]

# Square properties
square_size = 1200
square_x = (width - square_size) // 2
square_y = (height - square_size) // 2

# Animation settings
line_width = 1  # Width of the sweeping line (1 pixel)
sweep_speed = 5  # Speed of the sweep (1 pixel per frame)

def draw_square(animation_surface):
    pygame.draw.rect(animation_surface, (255, 255, 255, 0), (square_x, square_y, square_size, square_size), 1)  # Transparent fill

def sweep_lines(animation_surface):
    for x in range(square_x, square_x + square_size, sweep_speed):
        for color in colors:
            pygame.draw.line(animation_surface, color, (x, square_y), (x, square_y + square_size), line_width)
            screen.blit(background_image, (0, 0))  # Redraw background
            screen.blit(animation_surface, (0, 0))  # Apply the animation surface
            pygame.display.flip()  # Update the display
            pygame.time.wait(50)  # Wait a bit to see the color change before the next one

def main():
    running = True
    animation_surface = pygame.Surface((width, height), pygame.SRCALPHA)  # Use a surface with alpha channel

    screen.blit(background_image, (0, 0))  # Initial background blit
    draw_square(animation_surface)  # Draw the square outline
    sweep_lines(animation_surface)  # Sweep lines across the square

    # After the animation, keep the window open until the user decides to close it
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

    pygame.quit()
    sys.exit()

if __name__ == "__main__":
    main()
