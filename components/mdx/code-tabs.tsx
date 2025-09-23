"use client"

import { useState, type ReactNode } from "react"
import { cn } from "@/lib/utils"

interface CodeTabsProps {
  tabs: string[]
  children: ReactNode
}

export function CodeTabs({ tabs, children }: CodeTabsProps) {
  const [activeTab, setActiveTab] = useState(0)

  // Convert children to array and filter out non-code elements
  const childrenArray = Array.isArray(children) ? children : [children]
  const codeBlocks = childrenArray.filter(
    (child: any) => child?.props?.className?.includes("language-") || child?.type === "pre",
  )

  return (
    <div className="my-6 overflow-hidden rounded-lg border bg-card">
      <div className="flex border-b border-border bg-muted/50">
        {tabs.map((tab, index) => (
          <button
            key={tab}
            onClick={() => setActiveTab(index)}
            className={cn(
              "px-4 py-3 text-sm font-medium transition-colors relative",
              activeTab === index
                ? "bg-background text-foreground border-b-2 border-primary"
                : "text-muted-foreground hover:text-foreground hover:bg-muted/80",
            )}
          >
            {tab}
          </button>
        ))}
      </div>
      <div className="p-0">{codeBlocks[activeTab] || children}</div>
    </div>
  )
}
