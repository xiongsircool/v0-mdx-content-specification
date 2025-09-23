"use client"

import { useEffect, useState } from "react"

export function MetallicTitle() {
  const [timeOfDay, setTimeOfDay] = useState("morning")

  useEffect(() => {
    const updateTimeOfDay = () => {
      const hour = new Date().getHours()
      if (hour >= 5 && hour < 12) {
        setTimeOfDay("morning")
      } else if (hour >= 12 && hour < 17) {
        setTimeOfDay("afternoon")
      } else if (hour >= 17 && hour < 21) {
        setTimeOfDay("evening")
      } else {
        setTimeOfDay("night")
      }
    }

    updateTimeOfDay()
    const interval = setInterval(updateTimeOfDay, 60000) // Update every minute

    return () => clearInterval(interval)
  }, [])

  const getMetallicStyle = () => {
    const baseStyle = {
      backgroundSize: "200% 200%",
      WebkitBackgroundClip: "text",
      backgroundClip: "text",
      color: "#1a1a1a", // Fallback color for better visibility
      WebkitTextFillColor: "transparent", // This makes the gradient visible while keeping fallback
    }

    switch (timeOfDay) {
      case "morning":
        return {
          ...baseStyle,
          backgroundImage: "linear-gradient(135deg, #b8860b 0%, #ffd700 25%, #ffed4e 50%, #ffd700 75%, #b8860b 100%)",
        }
      case "afternoon":
        return {
          ...baseStyle,
          backgroundImage: "linear-gradient(135deg, #6b7280 0%, #9ca3af 25%, #d1d5db 50%, #9ca3af 75%, #6b7280 100%)",
        }
      case "evening":
        return {
          ...baseStyle,
          backgroundImage: "linear-gradient(135deg, #92400e 0%, #d97706 25%, #f59e0b 50%, #d97706 75%, #92400e 100%)",
        }
      case "night":
        return {
          ...baseStyle,
          backgroundImage: "linear-gradient(135deg, #374151 0%, #6b7280 25%, #9ca3af 50%, #6b7280 75%, #374151 100%)",
        }
      default:
        return baseStyle
    }
  }

  return (
    <h1
      className="text-6xl md:text-7xl font-bold mb-6 text-balance"
      style={{
        ...getMetallicStyle(),
        animation: "metallic-shine 3s ease-in-out infinite",
        textShadow: "0 2px 4px rgba(0,0,0,0.1)",
      }}
    >
      SBC 生物信息学俱乐部
    </h1>
  )
}
